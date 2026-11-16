from flask import Flask, render_template, request, jsonify, send_from_directory
from pathlib import Path
from pipeline import run_preprocessing

app = Flask(__name__)
ROOT = Path(__file__).resolve().parent
OUTPUTS = ROOT / "outputs"
OUTPUTS.mkdir(exist_ok=True)

@app.get("/")
def index():
    return render_template("index.html")

@app.post("/api/process")
def process():
    data = request.get_json(silent=True) or {}
    pdb_id = str(data.get("pdb_id", "")).strip().upper()
    chain_id = str(data.get("chain_id", "")).strip()
    try:
        if len(pdb_id) != 4 or not pdb_id.isalnum():
            return jsonify({"error": "Enter a valid four-character PDB ID."}), 400
        if not chain_id:
            return jsonify({"error": "Enter a chain ID."}), 400
        result = run_preprocessing(pdb_id, chain_id, OUTPUTS)
        return jsonify(result)
    except Exception as exc:
        app.logger.exception("Preprocessing failed")
        return jsonify({"error": str(exc)}), 500

@app.get("/download/<job>/<filename>")
def download(job, filename):
    folder = (OUTPUTS / job).resolve()
    if OUTPUTS.resolve() not in folder.parents:
        return "Invalid path", 400
    return send_from_directory(folder, filename, as_attachment=True)

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=True)
