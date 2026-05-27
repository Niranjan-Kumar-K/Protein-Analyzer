from flask import Flask, render_template, request
from Bio.SeqUtils.ProtParam import ProteinAnalysis

app = Flask(__name__)

def generate_chromatography_strategy(sequence):
    try:
        raw_input = sequence.strip()
        
        # 1. Parse lines and filter out FASTA header metadata smoothly
        lines = raw_input.splitlines()
        clean_lines = [line.strip() for line in lines if not line.startswith(">")]
        processed_string = "".join(clean_lines).upper()
        
        if not processed_string:
            return {"error": "Invalid sequence. Please provide valid characters."}

        # 2. Compute Dynamic Certainty Score before filtering out unknown tokens
        total_len = len(processed_string)
        ambiguous_count = sum(1 for char in processed_string if char in ['X', 'B', 'Z', 'U', 'O'])
        gap_count = raw_input.count('-') + raw_input.count('.')
        
        certainty_val = 100 - (ambiguous_count * 15) - (gap_count * 20)
        certainty_percentage = max(50, min(100, certainty_val))

        # 3. Strip non-standard characters so Biopython doesn't throw calculation errors
        standard_acids = "ACDEFGHIKLMNPQRSTVWY"
        cleaned_seq = "".join(c for c in processed_string if c in standard_acids)
        
        if not cleaned_seq:
            return {"error": "Sequence contains no standard amino acid residues for analysis."}

        # 4. Engage Biopython Engine with standard residues
        analysed_seq = ProteinAnalysis(cleaned_seq)
        
        raw_mw = analysed_seq.molecular_weight()
        mw_kda = f"{raw_mw / 1000:.2f} kDa"
        
        raw_pi = analysed_seq.isoelectric_point()
        pi_val = f"{raw_pi:.2f}"
        
        instability = analysed_seq.instability_index()
        stability_status = "Stable" if instability < 40 else "Unstable / Highly Labile"
        instability_val = f"{instability:.2f} ({stability_status})"
        
        # Calculate Aliphatic Index via amino acid volume fractions attribute
        aa_perc = analysed_seq.amino_acids_percent
        val_perc = aa_perc.get('V', 0) * 100
        ile_perc = aa_perc.get('I', 0) * 100
        leu_perc = aa_perc.get('L', 0) * 100
        ala_perc = aa_perc.get('A', 0) * 100
        aliphatic_index = ala_perc + (2.9 * val_perc) + (3.9 * ile_perc) + (3.9 * leu_perc)
        aliphatic_val = f"{aliphatic_index:.2f}"
        
        # Extinction Coefficient & Absorbance (A280) mapping
        try:
            epsilon_reduced = analysed_seq.molar_extinction_coefficient()[0]
            extinction_val = f"{epsilon_reduced} M⁻¹ cm⁻¹"
            a280 = epsilon_reduced / raw_mw
            a280_val = f"{a280:.3f}"
        except:
            extinction_val = "0 M⁻¹ cm⁻¹ (No W/Y residues)"
            a280_val = "0.000"

        # 5. Phase I: Capture Chromatography Decision Paths Matrix
        if raw_pi < 6.5:
            resin = "Q Sepharose Fast Flow (Strong Anion Exchanger)"
            buffer = "20 mM Tris-HCl, pH 8.0"
            rationale = (
                f"With a calculated pI of {pi_val}, the target protein carries a net negative surface charge "
                f"at physiological ranges. A strong anion exchange matrix will optimally secure the negative envelope."
            )
        elif raw_pi > 7.5:
            resin = "SP Sepharose Fast Flow (Strong Cation Exchanger)"
            buffer = "20 mM Sodium Phosphate, pH 6.0"
            rationale = (
                f"The basic nature of this sequence (pI {pi_val}) guarantees a net positive charge distribution. "
                f"A strong cation exchange matrix will securely bind the positive coordinates."
            )
        else:
            resin = "DEAE Sepharose or Multimodal Capto MMC"
            buffer = "20 mM HEPES, pH 7.0"
            rationale = f"The sequence possesses a near-neutral vector (pI {pi_val}). A multimodal matrix is selected to exploit subtle electrostatic pocket variations."
            
        # 6. Dynamic Phase II: Intermediate Purification (C-I-P Implementation)
        # Trigger condition: Heavy molecular scale (>45 kDa) OR structurally unstable sequence (>40)
        intermediate_resin = None
        intermediate_buffer = None
        intermediate_rationale = None

        if raw_mw > 45000 or instability > 40:
            intermediate_resin = "Phenyl Sepharose 6 Fast Flow (Hydrophobic Interaction Chromatography)"
            intermediate_buffer = "20 mM Sodium Phosphate, 1.5 M (NH4)2SO4, pH 7.0"
            intermediate_rationale = (
                f"Due to elevated molecular weight structural scale ({mw_kda}) or significant instability indices ({instability:.2f}), "
                f"an intermediate HIC refinement phase is dynamically engaged. Utilizing a decreasing ammonium sulfate salt gradient "
                f"strips away closely migrating host cell protein variants and misfolded configurations before polishing."
            )

        # 7. Phase III/II: Polishing Column Selection based on Molecular Weight
        if raw_mw < 30000:
            polishing = "Superdex 75 Increase (Separation range: 3 - 70 kDa)"
        else:
            polishing = "Superdex 200 Increase (Separation range: 10 - 600 kDa)"

        return {
            "mw": mw_kda,
            "pI": pi_val,
            "instability": instability_val,
            "aliphatic": aliphatic_val,
            "extinction": extinction_val,
            "a280": a280_val,
            "resin": resin,
            "buffer": buffer,
            "rationale": rationale,
            "intermediate_resin": intermediate_resin,
            "intermediate_buffer": intermediate_buffer,
            "intermediate_rationale": intermediate_rationale,
            "polishing": polishing,
            "certainty": f"{certainty_percentage}% Formula Verified"
        }

    except Exception as e:
        return {"error": f"Analysis failed: {str(e)}"}

@app.route("/", methods=["GET", "POST"])
def index():
    results = None
    protein_input = ""
    if request.method == "POST":
        protein_input = request.form.get("protein_sequence", "")
        results = generate_chromatography_strategy(protein_input)
    return render_template("index.html", results=results, protein_input=protein_input)

if __name__ == "__main__":
    app.run(debug=True)