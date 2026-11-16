const runBtn = document.getElementById("run");
runBtn.addEventListener("click", async () => {
  const pdb_id = document.getElementById("pdb").value.trim().toUpperCase();
  const chain_id = document.getElementById("chain").value.trim();
  const progress = document.getElementById("progress");
  const results = document.getElementById("results");
  results.classList.add("hidden");
  progress.classList.remove("hidden");
  runBtn.disabled = true;
  document.getElementById("progressText").textContent = "Downloading structure, extracting chain and running repair/validation…";
  try {
    const response = await fetch("/api/process", {
      method:"POST", headers:{"Content-Type":"application/json"},
      body:JSON.stringify({pdb_id, chain_id})
    });
    const data = await response.json();
    if (!response.ok) throw new Error(data.error || "Processing failed");
    render(data);
    progress.classList.add("hidden");
    results.classList.remove("hidden");
    results.scrollIntoView({behavior:"smooth"});
  } catch (err) {
    document.getElementById("progressText").textContent = "Error: " + err.message;
    document.getElementById("state").textContent = "Failed";
  } finally { runBtn.disabled = false; }
});

function render(data) {
  const checks = data.checks || [];
  const bad = checks.some(x => x.status === "FAIL");
  const warn = checks.some(x => x.status === "WARNING") || (data.warnings || []).length;
  document.getElementById("overall").textContent = bad ? "Failed" : warn ? "Completed with warnings" : "Validated";
  document.getElementById("summary").innerHTML = `
    <div class="metric"><strong>${data.output_structure.residues}</strong><span>Protein residues</span></div>
    <div class="metric"><strong>${data.output_structure.atom_records}</strong><span>Output atom records</span></div>
    <div class="metric"><strong>${checks.filter(x=>x.status==="PASS").length}/${checks.length}</strong><span>Checks passed</span></div>`;
  document.getElementById("checks").innerHTML = checks.map(x => `
    <div class="check"><i class="dot ${x.status.toLowerCase()}"></i><div><b>${escapeHtml(x.name)} · ${x.status}</b><p>${escapeHtml(x.detail)}</p></div></div>`).join("");
  document.getElementById("files").innerHTML = (data.files || []).map(f =>
    `<div class="file"><span>${escapeHtml(f.name)} <small>(${formatBytes(f.bytes)})</small></span><a href="${f.url}">Download ↓</a></div>`).join("");
  const warnings = data.warnings || [];
  const warningBox = document.getElementById("warnings");
  if (warnings.length) { warningBox.innerHTML = "<b>Scientific review required</b><br>" + warnings.map(escapeHtml).join("<br>"); warningBox.classList.remove("hidden"); }
  else warningBox.classList.add("hidden");
}
function escapeHtml(s){return String(s).replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));}
function formatBytes(n){return n<1024?n+" B":(n/1024).toFixed(1)+" KB";}
