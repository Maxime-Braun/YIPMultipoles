import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";


window.THREE = THREE;
window.OrbitControls = OrbitControls;

        // Modal Logic
        const helpModal = document.getElementById("helpModal");
        const helpBtn = document.getElementById("helpLink");
        const helpClose = document.getElementById("helpModalClose");

        const cellParamsModal = document.getElementById("cellParamsModal");
        const cellParamsClose = document.getElementById("cellParamsClose");
        const cellParamsApplyBtn = document.getElementById("applyCellParamsBtn");
        const cellParamsCancelBtn = document.getElementById("cancelCellParamsBtn");
        const cellParamsMsg = document.getElementById("cellParamsMsg");
        const cellParamsError = document.getElementById("cellParamsError");
        const importMcifBtn = document.getElementById("importMcifBtn");
        const mcifFileInput = document.getElementById("mcifFileInput");
        const cellParamInputs = {
          a: document.getElementById("cell_a"),
          b: document.getElementById("cell_b"),
          c: document.getElementById("cell_c"),
          alpha: document.getElementById("cell_alpha"),
          beta: document.getElementById("cell_beta"),
          gamma: document.getElementById("cell_gamma")
        };

        let CURRENT_CELL_PARAM_SYSTEM = null;
        const CUSTOM_CELL_PARAMS = {
          triclinic: { a: 1, b: 1, c: 1, alpha: 90, beta: 90, gamma: 90 },
          monoclinic: { a: 1.2, b: 1.0, c: 1.3, alpha: 90, beta: 110, gamma: 90 }
        };
  const CUSTOM_CELL_ENABLED = { triclinic: false, monoclinic: false };
  let CUSTOM_CELL_GLOBAL = null;
  let CUSTOM_CELL_GLOBAL_ENABLED = false;
  let MCIF_ATOM_VARS = null;
  let MCIF_ATOM_SITES = null;
  let SUPPRESS_CELL_PARAM_MODAL_ONCE = false;

        function parseMcifValue(raw) {
          if (!raw) return null;
          const cleaned = String(raw).replace(/[\"']/g, "").trim();
          const val = parseFloat(cleaned);
          return Number.isFinite(val) ? val : null;
        }

        function parseMcifString(raw) {
          if (!raw) return null;
          return String(raw).replace(/[\"']/g, "").trim();
        }

        function tokenizeMcifRow(line) {
          if (!line) return [];
          const tokens = line.match(/'(?:[^']*)'|"(?:[^"]*)"|\S+/g);
          return tokens ? tokens.map(t => t.trim()) : [];
        }

        function extractMcifRaw(text, tag) {
          const re = new RegExp(`${tag}\\s+([^\\r\\n]+)`, "i");
          const match = text.match(re);
          return match ? match[1].trim() : null;
        }

        function extractMcifTag(text, tag) {
          const raw = extractMcifRaw(text, tag);
          return raw ? parseMcifValue(raw) : null;
        }

        function parseMcifMagneticGroup(text) {
          const numRaw = extractMcifRaw(text, "_space_group_magn.number_BNS") ||
            extractMcifRaw(text, "_space_group_magn.number_bns") ||
            extractMcifRaw(text, "_magnetic_space_group.number_BNS");
          if (numRaw) {
            const match = numRaw.match(/[\d.]+/);
            if (match) return match[0];
          }

          const nameRaw = extractMcifRaw(text, "_space_group_magn.name_BNS") ||
            extractMcifRaw(text, "_space_group_magn.name_bns") ||
            extractMcifRaw(text, "_magnetic_space_group.name_BNS");
          if (nameRaw) {
            const bnsMatch = nameRaw.match(/BNS:\s*([\d.]+)/i);
            if (bnsMatch) return bnsMatch[1];
            const numMatch = nameRaw.match(/[\d.]+/);
            if (numMatch) return numMatch[0];
          }

          const fallback = text.match(/BNS:\s*([\d.]+)/i);
          if (fallback) return fallback[1];
          return null;
        }

        function parseMcifAtomFract(text) {
          const lines = text.split(/\r?\n/);
          let i = 0;
          while (i < lines.length) {
            const line = lines[i].trim();
            if (line.toLowerCase() === "loop_") {
              const headers = [];
              i += 1;
              while (i < lines.length) {
                const hdr = lines[i].trim();
                if (!hdr || hdr.startsWith("#")) {
                  i += 1;
                  continue;
                }
                if (!hdr.startsWith("_")) break;
                headers.push(hdr);
                i += 1;
              }
              const xIdx = headers.findIndex(h => h.toLowerCase() === "_atom_site_fract_x");
              const yIdx = headers.findIndex(h => h.toLowerCase() === "_atom_site_fract_y");
              const zIdx = headers.findIndex(h => h.toLowerCase() === "_atom_site_fract_z");
              if (xIdx >= 0 && yIdx >= 0 && zIdx >= 0) {
                while (i < lines.length) {
                  const dataLine = lines[i].trim();
                  if (!dataLine || dataLine.startsWith("#")) {
                    i += 1;
                    continue;
                  }
                  if (dataLine.startsWith("_") || dataLine.toLowerCase() === "loop_" || dataLine.startsWith("data_")) {
                    break;
                  }
                  const parts = tokenizeMcifRow(dataLine);
                  if (parts.length >= headers.length) {
                    const x = parseMcifValue(parts[xIdx]);
                    const y = parseMcifValue(parts[yIdx]);
                    const z = parseMcifValue(parts[zIdx]);
                    if ([x, y, z].every(v => v !== null)) {
                      return { x, y, z };
                    }
                  }
                  i += 1;
                }
              }
            }
            i += 1;
          }
          return null;
        }

        function parseMcifAtomSites(text) {
          const lines = text.split(/\r?\n/);
          let i = 0;
          while (i < lines.length) {
            const line = lines[i].trim();
            if (line.toLowerCase() === "loop_") {
              const headers = [];
              i += 1;
              while (i < lines.length) {
                const hdr = lines[i].trim();
                if (!hdr || hdr.startsWith("#")) {
                  i += 1;
                  continue;
                }
                if (!hdr.startsWith("_")) break;
                headers.push(hdr);
                i += 1;
              }

              const lowerHeaders = headers.map(h => h.toLowerCase());
              const labelIdx = lowerHeaders.findIndex(h =>
                h === "_atom_site_label" ||
                h === "_atom_site_type_symbol" ||
                h === "_atom_site_symbol"
              );
              const wyckoffIdx = lowerHeaders.findIndex(h =>
                h === "_atom_site_wyckoff_symbol" ||
                h === "_atom_site_wyckoff_letter" ||
                h === "_atom_site_wyckoff"
              );

              if (labelIdx >= 0 && wyckoffIdx >= 0) {
                const sites = [];
                while (i < lines.length) {
                  const dataLine = lines[i].trim();
                  if (!dataLine || dataLine.startsWith("#")) {
                    i += 1;
                    continue;
                  }
                  if (dataLine.startsWith("_") || dataLine.toLowerCase() === "loop_" || dataLine.startsWith("data_")) {
                    break;
                  }
                  const parts = tokenizeMcifRow(dataLine);
                  if (parts.length >= headers.length) {
                    const label = parseMcifString(parts[labelIdx]);
                    const wyckoff = parseMcifString(parts[wyckoffIdx]);
                    if (label && wyckoff) {
                      sites.push({ label, wyckoff });
                    }
                  }
                  i += 1;
                }

                if (sites.length > 0) return sites;
              }
            }
            i += 1;
          }
          return [];
        }

        function renderMcifAtomHeader(atomSites) {
          const header = document.getElementById("mcifAtomHeader");
          if (!header) return;

          if (!atomSites || atomSites.length === 0) {
            header.innerHTML = "";
            header.style.display = "none";
            return;
          }

          const items = atomSites
            .map(site => `${site.label}: ${site.wyckoff}`)
            .join(", ");
          header.textContent = `Atoms (Wyckoff): ${items}`;
          header.style.display = "block";
        }

        function parseMcifCellParams(text) {
          const params = {
            a: extractMcifTag(text, "_cell_length_a"),
            b: extractMcifTag(text, "_cell_length_b"),
            c: extractMcifTag(text, "_cell_length_c"),
            alpha: extractMcifTag(text, "_cell_angle_alpha"),
            beta: extractMcifTag(text, "_cell_angle_beta"),
            gamma: extractMcifTag(text, "_cell_angle_gamma")
          };
          if (Object.values(params).some(v => v === null)) return null;
          return params;
        }

        function showCellParamsError(msg) {
          if (!cellParamsError) return;
          cellParamsError.textContent = msg;
          cellParamsError.style.display = msg ? "block" : "none";
        }

        function readCellParamsFromInputs() {
          return {
            a: parseFloat(cellParamInputs.a?.value),
            b: parseFloat(cellParamInputs.b?.value),
            c: parseFloat(cellParamInputs.c?.value),
            alpha: parseFloat(cellParamInputs.alpha?.value),
            beta: parseFloat(cellParamInputs.beta?.value),
            gamma: parseFloat(cellParamInputs.gamma?.value)
          };
        }

        function validateCellParams(params) {
          const vals = Object.values(params);
          if (vals.some(v => !isFinite(v))) {
            return "Please enter numeric values for all cell parameters.";
          }
          if (params.a <= 0 || params.b <= 0 || params.c <= 0) {
            return "Cell lengths must be positive.";
          }
          if (params.alpha <= 0 || params.beta <= 0 || params.gamma <= 0 ||
              params.alpha >= 180 || params.beta >= 180 || params.gamma >= 180) {
            return "Angles must be between 0° and 180° (exclusive).";
          }
          const gammaRad = params.gamma * Math.PI / 180;
          if (Math.abs(Math.sin(gammaRad)) < 1e-6) {
            return "Gamma angle cannot be 0° or 180° for a valid basis.";
          }
          return null;
        }

        function setCellParamInputs(system) {
          const base = CUSTOM_CELL_PARAMS[system] || { a: 1, b: 1, c: 1, alpha: 90, beta: 90, gamma: 90 };
          if (cellParamInputs.a) cellParamInputs.a.value = base.a;
          if (cellParamInputs.b) cellParamInputs.b.value = base.b;
          if (cellParamInputs.c) cellParamInputs.c.value = base.c;
          if (cellParamInputs.alpha) cellParamInputs.alpha.value = base.alpha;
          if (cellParamInputs.beta) cellParamInputs.beta.value = base.beta;
          if (cellParamInputs.gamma) cellParamInputs.gamma.value = base.gamma;
        }

        function setCellParamsModalMode(showInputs) {
          const grid = document.getElementById("cellParamsGrid");
          const actions = document.getElementById("cellParamsActions");
          if (grid) grid.style.display = showInputs ? "grid" : "none";
          if (actions) actions.style.display = showInputs ? "flex" : "none";
        }

        function openCellParamsModal(system) {
          if (!cellParamsModal) return;
          CURRENT_CELL_PARAM_SYSTEM = system;
          setCellParamsModalMode(true);
          if (cellParamsMsg) {
            if (system === "monoclinic") {
              cellParamsMsg.textContent =
                "Monoclinic system detected. Enter lattice parameters (α, γ typically 90°).";
            } else if (system === "hexagonal" || system === "trigonal") {
              cellParamsMsg.textContent =
                "Hexagonal setting required. Enter hexagonal cell parameters (a=b, α=β=90°, γ=120°).";
            } else {
              cellParamsMsg.textContent =
                "Triclinic system detected. Enter full lattice parameters.";
            }
          }
          setCellParamInputs(system);
          showCellParamsError("");
          cellParamsModal.style.display = "block";
        }

        function openCellParamsErrorModal(message) {
          if (!cellParamsModal) return;
          CURRENT_CELL_PARAM_SYSTEM = null;
          setCellParamsModalMode(false);
          if (cellParamsMsg) {
            cellParamsMsg.textContent = message;
          }
          showCellParamsError("");
          cellParamsModal.style.display = "block";
        }

        function closeCellParamsModal() {
          if (!cellParamsModal) return;
          cellParamsModal.style.display = "none";
          showCellParamsError("");
          setCellParamsModalMode(true);
        }

        function maybeOpenCellParamsModal(system) {
          if (system !== "triclinic" && system !== "monoclinic") return;
          if (SUPPRESS_CELL_PARAM_MODAL_ONCE) {
            SUPPRESS_CELL_PARAM_MODAL_ONCE = false;
            return;
          }
          if (cellParamsModal && cellParamsModal.style.display === "block") return;
          openCellParamsModal(system);
        }

        function applyCellParams(system, params) {
          if (system && system !== "triclinic" && system !== "monoclinic") {
            applyGlobalCellParams(params);
            closeCellParamsModal();
            return;
          }
          CUSTOM_CELL_PARAMS[system] = { ...params };
          CUSTOM_CELL_ENABLED[system] = true;
          if (currentGroup && currentGroup.name) {
            SUPPRESS_CELL_PARAM_MODAL_ONCE = true;
            setUnitCellFromGroupName(currentGroup.name);
          }
          closeCellParamsModal();
        }
        function applyGlobalCellParams(params) {
          CUSTOM_CELL_GLOBAL = { ...params };
          CUSTOM_CELL_GLOBAL_ENABLED = true;
          if (currentGroup && currentGroup.name) {
            SUPPRESS_CELL_PARAM_MODAL_ONCE = true;
            setUnitCellFromGroupName(currentGroup.name);
          }
        }

        helpBtn.onclick = function(e) {
            e.preventDefault();
            helpModal.style.display = "block";
        }
        helpClose.onclick = function() {
            helpModal.style.display = "none";
        }

        if (importMcifBtn && mcifFileInput) {
          importMcifBtn.onclick = () => mcifFileInput.click();
          mcifFileInput.onchange = async () => {
            const file = mcifFileInput.files && mcifFileInput.files[0];
            if (!file) return;
            const text = await file.text();
            const params = parseMcifCellParams(text);
            if (!params) {
              const status = document.getElementById("status");
              if (status) status.innerText = "Could not read cell parameters from the mcif file.";
              renderMcifAtomHeader(null);
              mcifFileInput.value = "";
              return;
            }

            const status = document.getElementById("status");
            const groupInput = document.getElementById("groupNumberInput");
            const mcifGroup = parseMcifMagneticGroup(text);
            const mcifAtom = parseMcifAtomFract(text);
            const mcifAtomSites = parseMcifAtomSites(text);
            MCIF_ATOM_VARS = mcifAtom ? { ...mcifAtom } : null;
            MCIF_ATOM_SITES = mcifAtomSites ? [...mcifAtomSites] : null;
            renderMcifAtomHeader(MCIF_ATOM_SITES);

            if (mcifGroup && groupInput) {
              groupInput.value = mcifGroup;
              if (typeof groupInput.onchange === "function") {
                await groupInput.onchange({ target: groupInput });
              } else {
                groupInput.dispatchEvent(new Event("change"));
              }
            }

            const system = currentGroup && currentGroup.name
              ? crystalSystemFromSGNumber(parseBNSNumberMajor(currentGroup.name))
              : null;

            if (system === "hexagonal" || system === "trigonal") {
              const near = (a, b, eps = 1e-2) => Math.abs(a - b) < eps;
              const isHexCell =
                near(params.a, params.b) &&
                near(params.alpha, 90) &&
                near(params.beta, 90) &&
                near(params.gamma, 120);
              if (!isHexCell) {
                if (status) {
                  status.innerText =
                    "Group selected, but cell parameters were rejected. Please provide hexagonal cell parameters.";
                }
                openCellParamsErrorModal(
                  "Hexagonal setting required. This mcif uses a rhombohedral/non-hexagonal cell. Please provide hexagonal cell parameters (a=b, α=β=90°, γ=120°)."
                );
                renderMcifAtomHeader(null);
                mcifFileInput.value = "";
                return;
              }
            }

            CUSTOM_CELL_PARAMS.triclinic = { ...params };
            CUSTOM_CELL_PARAMS.monoclinic = { ...params };
            applyGlobalCellParams(params);

            if (system === "triclinic" || system === "monoclinic") {
              setCellParamInputs(system);
              if (cellParamInputs.a) cellParamInputs.a.value = params.a;
              if (cellParamInputs.b) cellParamInputs.b.value = params.b;
              if (cellParamInputs.c) cellParamInputs.c.value = params.c;
              if (cellParamInputs.alpha) cellParamInputs.alpha.value = params.alpha;
              if (cellParamInputs.beta) cellParamInputs.beta.value = params.beta;
              if (cellParamInputs.gamma) cellParamInputs.gamma.value = params.gamma;
              const err = validateCellParams(params);
              if (err) {
                showCellParamsError(err);
                openCellParamsModal(system);
              } else {
                applyCellParams(system, params);
                if (status) {
                  status.innerText = mcifGroup
                    ? `Cell parameters loaded from mcif. Group ${mcifGroup} selected.`
                    : "Cell parameters loaded from mcif.";
                }
              }
            } else {
              if (status) {
                status.innerText = mcifGroup
                  ? `Group ${mcifGroup} selected. Cell parameters loaded from mcif.`
                  : "Cell parameters loaded from mcif.";
              }
            }
            mcifFileInput.value = "";
          };
        }

        if (cellParamsClose) {
          cellParamsClose.onclick = function() {
            closeCellParamsModal();
          }
        }
        if (cellParamsCancelBtn) {
          cellParamsCancelBtn.onclick = function() {
            if (CURRENT_CELL_PARAM_SYSTEM) {
              CUSTOM_CELL_ENABLED[CURRENT_CELL_PARAM_SYSTEM] = false;
              if (currentGroup && currentGroup.name) {
                SUPPRESS_CELL_PARAM_MODAL_ONCE = true;
                setUnitCellFromGroupName(currentGroup.name);
              }
            }
            closeCellParamsModal();
          }
        }
        if (cellParamsApplyBtn) {
          cellParamsApplyBtn.onclick = function() {
            if (!CURRENT_CELL_PARAM_SYSTEM) return;
            const params = readCellParamsFromInputs();
            const err = validateCellParams(params);
            if (err) {
              showCellParamsError(err);
              return;
            }
            applyCellParams(CURRENT_CELL_PARAM_SYSTEM, params);
          }
        }

        window.onclick = function(event) {
            if (event.target == helpModal) {
                helpModal.style.display = "none";
            }
      if (event.target == cellParamsModal) {
        closeCellParamsModal();
            }
        }
// ==================== Tensor Ylm state ====================
// rank → { harmonics:[[l,m],...], weights:number[] }
        let CELL_BASIS = { a:[1,0,0], b:[0,1,0], c:[0,0,1] };

        let db = [];
        let currentGroup = null;
        let lastResults = []; // Store results for toggling view
        window.lastResults = lastResults;
        let wyckoffRandomVars = null; // {x,y,z} used for the current Wyckoff orbit
  const MIN_ORBIT_DISTANCE = 0.12; // in ideal-cartesian units (cell length ~ 1)
  const MAX_ORBIT_TRIES = 60;
        // -------- 3D Viewer --------
        let scene3d, camera3d, renderer3d, controls3d;
        let orbitGroup3d = null; // holds the spheres for current orbit
        let axesHelper3d = null;
        // ==================== Atom hover tooltip ====================
        let atomTooltipEl = null;
        let atomRaycaster = null;
        let atomMouseNdc = null;
        let __atomHoverBound = false;


    // Isosurface meshes group (rank visualization)
    let isoGroup3d = null;
    let crystalAxesGroup3d = null;
    let SHOW_CRYSTAL_AXES = true;

    // Legacy hook: older tensor_visualizer UI had Ylm sliders.
    // In this app we expose the site tensor using the SAME Y(l,m) slider UI as tensor_visualizer (built from the site-symmetry basis).
    // Keep a no-op to avoid runtime errors if compute() calls it.
    // If true: backend already provides per-site tensors/multipoles in global frame.
// Rendering must NOT apply extra symmetry sign (detR*tr) per atom.
    let TENSORS_ALREADY_PER_SITE = true;
    let multipoleMode = 'magnetic'; // 'magnetic' or 'electric'
// Global registry of canonical variable labels (persistent across atoms in current position)
const VARIABLE_REGISTRY = {};

        
// Load Data (moved to backend API)
const API_BASE = (import.meta.env.VITE_API_BASE || (window.location.origin + "/api"));
        async function fetchJSON(url, options) {
            const res = await fetch(url, options);
            if (!res.ok) {
                const txt = await res.text();
                throw new Error(`HTTP ${res.status}: ${txt.substring(0,200)}`);
            }
            return await res.json();
        }

        async function loadGroupsFromBackend() {
            const groups = await fetchJSON(`${API_BASE}/groups`);
            db = groups; // [{index, name}]
            const sel = document.getElementById('groupSelect');
            sel.innerHTML = '';
      const placeholder = document.createElement('option');
      placeholder.value = '';
      placeholder.text = 'Select MSG...';
      sel.appendChild(placeholder);
            db.forEach(g => {
                const opt = document.createElement('option');
                opt.value = g.index;
                opt.text = g.name;
                sel.appendChild(opt);
            });

            const updateInputFromSelect = () => {
                const idx = parseInt(sel.value);
                const g = db.find(x => x.index === idx);
                if (!g) return;
                const match = (g.name || "").match(/BNS:\s*([\d\.]+)/);
                document.getElementById('groupNumberInput').value = match ? match[1] : '';
            };

            sel.onchange = async () => {
            updateInputFromSelect();

            resetWyckoffSelectionAndClearViewer();
            document.getElementById("results").innerHTML = "";
            document.getElementById("status").innerText = "Select a Wyckoff position…";

            await populateWyckoff(sel.value);

            const wSel = document.getElementById("wyckoffSelect");
            if (wSel) {
                wSel.value = "";
                wSel.selectedIndex = 0;
            }
            };

            document.getElementById('groupNumberInput').onchange = async (e) => {
                const val = (e.target.value || '').trim();
                if (!val) return;

                let idx = db.findIndex(g => {
                    const match = (g.name || "").match(/BNS:\s*([\d\.]+)/);
                    return match && match[1] === val;
                });

                if (idx === -1) {
                    idx = db.findIndex(g => {
                        const match = (g.name || "").match(/BNS:\s*([\d\.]+)/);
                        return match && match[1] && match[1].startsWith(val);
                    });
                }

               if (idx !== -1) {
                    // BUG 4 FIX: group changed => clear everything
                    resetWyckoffSelectionAndClearViewer();
                    document.getElementById("results").innerHTML = "";
                    document.getElementById("status").innerText = "Select a Wyckoff position…";

                    sel.value = db[idx].index;
                    updateInputFromSelect();
                    await populateWyckoff(sel.value);

                    // force placeholder (so nothing appears until user selects)
                    const wSel = document.getElementById("wyckoffSelect");
                    if (wSel) { wSel.value = ""; wSel.selectedIndex = 0; }
                }

            };

      sel.value = '';
        }

        async function loadGroupDetails(idx) {
            return await fetchJSON(`${API_BASE}/group/${idx}`);
        }
// ==================== Spin channel (tensor_visualizer convention) ====================
// 0->Mx, 1->My, 2->Mz
let ACTIVE_SPIN = 0;

function ensureSpinButtons() {
  const box = document.getElementById('tensor-sliders');
  if (!box) return;

  if (document.getElementById('site-spin-row')) return;

  const row = document.createElement('div');
  row.id = 'site-spin-row';
  row.style.display = 'flex';
  row.style.gap = '6px';
  row.style.alignItems = 'center';
  row.style.margin = '6px 0 10px 0';

  const lab = document.createElement('div');
  lab.textContent = 'Channel:';
  lab.style.fontWeight = '700';
  row.appendChild(lab);

  const names = ['Mx','My','Mz'];

  function refreshButtons() {
    Array.from(row.querySelectorAll('button')).forEach((b, i) => {
      b.style.background = (i === ACTIVE_SPIN) ? '#e8eefc' : '#fff';
    });
  }

  names.forEach((name, i) => {
    const b = document.createElement('button');
    b.type = 'button';
    b.textContent = name;
    b.style.padding = '4px 10px';
    b.style.border = '1px solid #ccc';
    b.style.borderRadius = '6px';
    b.style.cursor = 'pointer';
    b.onclick = () => {
    ACTIVE_SPIN = i;
    refreshButtons();

    // keep current Wyckoff selection
    document.getElementById("status").innerText = "Spin changed.";

    // If we already have an orbit rendered, just re-render tensors/sliders for it
    if (_lastOrbitForTensor) {
        renderTensorIsosurfacesFromSiteTensor(_lastOrbitForTensor);
    } else {
        // If no orbit cached yet, try to render from current selection (if any)
        updateOrbitViewer();
    }
    };
    row.appendChild(b);
  });

  box.prepend(row);

  if (!document.getElementById('reset-sliders-btn')) {
    const resetBtn = document.createElement('button');
    resetBtn.id = 'reset-sliders-btn';
    resetBtn.type = 'button';
    resetBtn.textContent = 'Reset';
    resetBtn.style.padding = '4px 10px';
    resetBtn.style.border = '1px solid #ccc';
    resetBtn.style.borderRadius = '6px';
    resetBtn.style.cursor = 'pointer';
    resetBtn.style.marginLeft = '6px';
    resetBtn.onclick = () => {
      // set all sliders for the current rank/spin to zero
      const rank = getSelectedRank();
      if (rank) {
        if (GLOBAL_BASIS_COEFFS[rank]) {
          GLOBAL_BASIS_COEFFS[rank] = GLOBAL_BASIS_COEFFS[rank].map(() => 0);
        }
      }

      // update slider UI values
      document.querySelectorAll('#tensor-sliders input[type="range"]').forEach(sl => {
        sl.value = 0;
        const valEl = sl.parentElement?.querySelector('div[style*="width: 55px"]');
        if (valEl) valEl.textContent = '0.00';
      });

      if (rank) {
        const basis = getBasisForRank(rank);
        if (basis && basis.length) {
          scheduleIsoRender(rank, basis, ACTIVE_SPIN);
        } else {
          updateOrbitViewer();
        }
      } else {
        updateOrbitViewer();
      }
    };
    row.appendChild(resetBtn);
  }
  refreshButtons();
}

document.addEventListener('DOMContentLoaded', ensureSpinButtons);
function sliceSpinComponent(fullTensorFlat, fullRank, spinIndex) {
  // fullTensorFlat length = 3^fullRank
  // convention: first index is spin (0,1,2), remaining (fullRank-1) are spatial
  const spatialRank = fullRank - 1;
  const block = Math.pow(3, spatialRank);
  const offset = spinIndex * block;

  // return the spatial tensor slice of length 3^(fullRank-1)
  return fullTensorFlat.slice(offset, offset + block);
}
function clearCrystalAxes() {
  if (!scene3d) return;
  if (crystalAxesGroup3d) {
    scene3d.remove(crystalAxesGroup3d);
    crystalAxesGroup3d = null;
  }
}

function addCrystalAxes({ scale = 0.6, headLength = 0.12, headWidth = 0.08 } = {}) {
  if (!scene3d) return;

  clearCrystalAxes();
  if (!SHOW_CRYSTAL_AXES) return;

  crystalAxesGroup3d = new THREE.Group();

  const { a, b, c } = CELL_BASIS || { a:[1,0,0], b:[0,1,0], c:[0,0,1] };

  // origin is at the CENTER of the cell in your viewer convention (u-0.5, v-0.5, w-0.5) :contentReference[oaicite:2]{index=2}
  const origin = new THREE.Vector3(0, 0, 0);

  function arrowFromVec(vec3, colorHex, labelText) {
    const dir = vec3.clone();
    const len = dir.length();
    if (len < 1e-12) return;

    dir.normalize();

    const arrow = new THREE.ArrowHelper(
      dir,
      origin,
      len * scale,
      colorHex,
      headLength,
      headWidth
    );

    // add label near tip
    const tip = origin.clone().add(dir.clone().multiplyScalar(len * scale * 1.05));
    const label = makeTextSprite(labelText, colorHex);
    label.position.copy(tip);

    crystalAxesGroup3d.add(arrow);
    crystalAxesGroup3d.add(label);
  }

  arrowFromVec(new THREE.Vector3(a[0], a[1], a[2]), 0xff0000, "a");
  arrowFromVec(new THREE.Vector3(b[0], b[1], b[2]), 0x00aa00, "b");
  arrowFromVec(new THREE.Vector3(c[0], c[1], c[2]), 0x0000ff, "c");

  scene3d.add(crystalAxesGroup3d);
}

function toggleCrystalAxes() {
  SHOW_CRYSTAL_AXES = !SHOW_CRYSTAL_AXES;
  addCrystalAxes({ scale: 0.45 }); // redraw with new state
}

// Simple text sprite (no extra libs needed)
function makeTextSprite(text, colorHex) {
  const canvas = document.createElement("canvas");
  const ctx = canvas.getContext("2d");

  // crisp text
  const fontSize = 64;
  ctx.font = `bold ${fontSize}px monospace`;
  const padding = 18;
  const w = Math.ceil(ctx.measureText(text).width + padding * 2);
  const h = fontSize + padding * 2;

  canvas.width = w;
  canvas.height = h;

  ctx.font = `bold ${fontSize}px monospace`;
  ctx.textBaseline = "middle";
  ctx.textAlign = "center";

  // background

  // text
  const hex = "#" + (colorHex >>> 0).toString(16).padStart(6, "0");
  ctx.fillStyle = hex;
  ctx.fillText(text, w / 2, h / 2);

  const tex = new THREE.CanvasTexture(canvas);
  tex.minFilter = THREE.LinearFilter;
  tex.magFilter = THREE.LinearFilter;

  const mat = new THREE.SpriteMaterial({ map: tex, transparent: true, depthTest: false });
  const sprite = new THREE.Sprite(mat);

  // size in world units (tweak if you want bigger)
  const worldScale = 0.18;
  sprite.scale.set(worldScale * (w / h), worldScale, 1);

  return sprite;
}

function renderWyckoffOccupancyLine() {
  const boxWrap = document.getElementById("wyckoffOccBoxes");
  if (!boxWrap || !currentGroup || !Array.isArray(currentGroup.wyckoff)) return;

  boxWrap.innerHTML = "";
  // All Wyckoff positions (display only)
  currentGroup.wyckoff.forEach(w => {
    const item = document.createElement("span");
    item.className = "wyckoff-item";

    const cb = document.createElement("input");
    cb.type = "checkbox";
    cb.checked = false; // ALWAYS unchecked, no logic

    const txt = document.createElement("span");
    txt.textContent = w.label;

    item.appendChild(cb);
    item.appendChild(txt);
    boxWrap.appendChild(item);
  });
}

function fracToCartFromMatrix(u, v, w, M) {
  return [
    M[0][0]*u + M[0][1]*v + M[0][2]*w,
    M[1][0]*u + M[1][1]*v + M[1][2]*w,
    M[2][0]*u + M[2][1]*v + M[2][2]*w,
  ];
}

function minOrbitDistance(orbit, group) {
  if (!orbit || orbit.length < 2 || !group) return Infinity;
  const M = getIdealLatticeMatrix(group);
  let minDist = Infinity;
  for (let i = 0; i < orbit.length; i++) {
    const p1 = orbit[i].coord;
    for (let j = i + 1; j < orbit.length; j++) {
      const p2 = orbit[j].coord;
      const diff = [
        p2[0] - p1[0],
        p2[1] - p1[1],
        p2[2] - p1[2]
      ].map(d => d - Math.round(d)); // minimal image in fractional space

      const cart = fracToCartFromMatrix(diff[0], diff[1], diff[2], M);
      const dist = Math.sqrt(cart[0]*cart[0] + cart[1]*cart[1] + cart[2]*cart[2]);
      if (dist < minDist) minDist = dist;
    }
  }
  return minDist;
}

function generateWyckoffRandomVars(group, wyckoff) {
  if (MCIF_ATOM_VARS) return { ...MCIF_ATOM_VARS };
  if (!group || !wyckoff) return { x: Math.random(), y: Math.random(), z: Math.random() };

  let bestVars = null;
  let bestMin = -1;

  for (let i = 0; i < MAX_ORBIT_TRIES; i++) {
    const vars = { x: Math.random(), y: Math.random(), z: Math.random() };
    const orbit = getWyckoffOrbit(group, wyckoff, vars);
    const dmin = minOrbitDistance(orbit, group);
    if (dmin > bestMin) {
      bestMin = dmin;
      bestVars = vars;
    }
    if (dmin >= MIN_ORBIT_DISTANCE) return vars;
  }

  return bestVars || { x: Math.random(), y: Math.random(), z: Math.random() };
}
async function populateWyckoff(idx) {
            currentGroup = await loadGroupDetails(idx);
            setUnitCellFromGroupName(currentGroup.name);
            const sel = document.getElementById('wyckoffSelect');
            sel.innerHTML = '';
            
            // Add placeholder
            const placeholder = document.createElement('option');
            placeholder.value = "";
            placeholder.text = "Select Position...";
            sel.appendChild(placeholder);

            // Add option to compute using the whole space group (all symmetry operations)
            const wholeOpt = document.createElement('option');
            wholeOpt.value = "whole";
            wholeOpt.text = "System";
            sel.appendChild(wholeOpt);

            currentGroup.wyckoff.forEach((w, i) => {
                const opt = document.createElement('option');
                opt.value = i;
                opt.text = `${w.label} (${w.coord})`;
                sel.appendChild(opt);
            });
            
            // Clear previous results
            document.getElementById('results').innerHTML = '';
            document.getElementById('status').innerText = '';
            renderWyckoffOccupancyLine();
            // Automatically compute (will do nothing if placeholder is selected)
            
        }

    const FIXED_MAX_RANK = 5;

    const triggerCompute = async () => {
            const wSel = document.getElementById("wyckoffSelect")?.value;
            if (!wSel) {
            resetWyckoffSelectionAndClearViewer();
            return;
            }
            const wIdx = document.getElementById('wyckoffSelect').value;       // this clears #wyckoffBoxes   // <- ADD THIS LINE (rebuild it)

      const maxRank = FIXED_MAX_RANK;
      syncRankIsoWithMaxRank(maxRank);
            
            if (!currentGroup || wIdx === "") return;
            // Generate one random (x,y,z) per Wyckoff choice (not for "whole")
            if (wIdx !== 'whole') {
        wyckoffRandomVars = generateWyckoffRandomVars(currentGroup, currentGroup.wyckoff[wIdx]);
            } else {
                wyckoffRandomVars = null;
            }
            let wyckoff = null;
            if (wIdx !== 'whole') wyckoff = currentGroup.wyckoff[wIdx];

            try {
                lastResults = [];
                window.lastResults = lastResults;
                await compute(currentGroup, wyckoff, maxRank);
                window.lastResults = lastResults;
                if (wIdx !== "whole" && wyckoff) {
                    if (!wyckoffRandomVars) {
                        wyckoffRandomVars = { x: Math.random(), y: Math.random(), z: Math.random() };
                    }
                    const orbit = getWyckoffOrbit(currentGroup, wyckoff, wyckoffRandomVars);
                    renderOrbitSpheres(orbit);
                    } else {
                    clearViewerOrbit();
                    }

            } catch (e) {
                console.error(e);
                document.getElementById('status').innerText = "Error: " + e.message;
            }
        };

        document.getElementById('wyckoffSelect').onchange = () => {
        wyckoffRandomVars = null;   // new random for new wyckoff
        triggerCompute();           // your normal compute
        updateOrbitViewer();        // draw spheres
        };
        document.getElementById('voigtCheck').onchange = () => {
            const resultsDiv = document.getElementById('results');
            // Keep the site symmetry block (first child)
            const opsDiv = resultsDiv.querySelector('.rank-block'); 
            resultsDiv.innerHTML = '';
            if (opsDiv) resultsDiv.appendChild(opsDiv);
            
            lastResults.forEach(res => {
                displayResult(res.rank, res.basis, resultsDiv, res.group, res.wyckoff);
            });
        };

    // Header switch: show/hide all label definitions blocks
    const showLabelDefsSwitch = document.getElementById('showLabelDefsSwitch');
    if (showLabelDefsSwitch) {
      showLabelDefsSwitch.onchange = () => {
        const els = document.querySelectorAll('.label-defs');
        els.forEach(el => { el.style.display = showLabelDefsSwitch.checked ? 'block' : 'none'; });
      };
    }

    const includeSocSwitch = document.getElementById('includeSocSwitch');
    if (includeSocSwitch) {
      includeSocSwitch.onchange = async () => {
        const wSel = document.getElementById("wyckoffSelect")?.value;
        if (multipoleMode !== 'magnetic') {
          document.getElementById("status").innerText = "SOC toggle applies to Magnetic mode only.";
          return;
        }
        if (wSel) {
          document.getElementById("status").innerText = includeSocSwitch.checked
            ? "SOC enabled — recomputing…"
            : "SOC disabled — computing No-SOC basis…";
          resetTensorStateForNewPosition();
          await triggerCompute();
          updateOrbitViewer();
        }
      };
    }

    // Multipole type select menu
    const multipoleTypeSelect = document.getElementById('multipoleTypeSelect');
    if (multipoleTypeSelect) {
      if (multipoleTypeSelect.value !== multipoleMode) {
        multipoleTypeSelect.value = multipoleMode;
      }
      multipoleMode = multipoleTypeSelect.value || multipoleMode;
      multipoleTypeSelect.onchange = async () => {
        multipoleMode = multipoleTypeSelect.value;
        resetTensorStateForNewPosition();

        const wSel = document.getElementById("wyckoffSelect")?.value;
        if (wSel) {
          document.getElementById("status").innerText =
            "Mode changed — recomputing…";
          await triggerCompute();
          updateOrbitViewer();
        } else {
          resetTensorSliderMemory();
          clearIsoGroup();
          document.getElementById("status").innerText =
            "Mode changed — select a Wyckoff position…";
        }
      };
    }

        

        // --- Math & Physics Engine ---

        // --- Math & Physics Engine (moved to backend) ---

        async function compute(group, wyckoff, maxRank) {
            const status = document.getElementById('status');
            const resultsDiv = document.getElementById('results');
            resultsDiv.innerHTML = '';
            status.innerText = "Computing site symmetry...";

            // Yield to UI
            await new Promise(r => setTimeout(r, 10));

            const group_index = parseInt(document.getElementById('groupSelect').value);
            const wSel = document.getElementById('wyckoffSelect').value;

            const voigt = document.getElementById('voigtCheck').checked;
            const crystal_basis = document.getElementById('crystalBasisSwitch').checked;

            // Backend computes invariant basis; frontend keeps identical rendering.
            const payload = {
                group_index: group_index,
                wyckoff_index: wSel,
                max_rank: maxRank,
                mode: multipoleMode,
                voigt: voigt,
        crystal_basis: crystal_basis,
        include_soc: document.getElementById('includeSocSwitch') ? document.getElementById('includeSocSwitch').checked : true
            };

            const data = await fetchJSON(`${API_BASE}/compute`, {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify(payload)
            });

      const includeSoc = document.getElementById('includeSocSwitch') ? document.getElementById('includeSocSwitch').checked : true;
      const noSocForceForbidden = !includeSoc && multipoleMode === 'magnetic' &&
        (data.ranks || []).some(r => r.rank === 1 && (!r.basis || r.basis.length === 0));
      if (noSocForceForbidden) {
        clearIsoGroup();
        clearTensorSlidersUI();
      }

            // Display selected symmetry operations (identical HTML structure to original)
            const siteOps = data.site_ops || [];
            const opsDiv = document.createElement('div');
            opsDiv.className = 'rank-block';
            const opsTitle = (wSel !== 'whole') ? `Site Symmetry Operations (${siteOps.length})` : `Space Group Symmetry Operations (${siteOps.length})`;
            opsDiv.innerHTML = `<div class="rank-title">${opsTitle}</div>`;
            let opsHtml = '<div style="display: grid; grid-template-columns: repeat(auto-fill, minmax(250px, 1fr)); gap: 5px; padding: 10px; font-family: monospace;">';
            siteOps.forEach((op, i) => {
                let s = op.str || "Unknown";
                if (op.tr === -1) s += " (Time Rev.)";
                opsHtml += `<div>${i+1}. ${s}</div>`;
            });
            opsHtml += '</div>';
            opsDiv.innerHTML += opsHtml;
            resultsDiv.appendChild(opsDiv);

            // Fill ranks (basis vectors) and render using original displayResult()
            for (let rank = 1; rank <= maxRank; rank++) {
                const r = (data.ranks || []).find(x => x.rank === rank);
        const basis = (noSocForceForbidden || !r) ? [] : r.basis;
                lastResults.push({rank: rank, basis: basis, group: group, wyckoff: wyckoff});
                // --- Ylm slider pruning: union over basis vectors (tensor_visualizer logic)
                displayResult(rank, basis, resultsDiv, group, wyckoff);
            }

            // Rebuild Ylm sliders now that we know which harmonics are allowed by symmetry
    
            status.innerText = "Done.";
        }


// [Removed compute - moved to backend]


        
// [Removed computeMultipoles - moved to backend]


        function applyOpToTensor(currentV, rank, R, nextV) {
            // currentV is input (and current buffer). nextV is swap buffer.
            // We iterate rank times.
            
            const stride = [1];
            for(let i=1; i<=rank; i++) stride.push(stride[i-1]*3);
            
            for (let k = 0; k < rank; k++) {
                const axisStride = stride[rank - 1 - k];
                
                const blockSize = axisStride * 3;
                const numBlocks = currentV.length / blockSize;
                
                for (let bIdx = 0; bIdx < numBlocks; bIdx++) {
                    const blockStart = bIdx * blockSize;
                    for (let offset = 0; offset < axisStride; offset++) {
                        const idx0 = blockStart + offset;
                        const idx1 = idx0 + axisStride;
                        const idx2 = idx1 + axisStride;
                        
                        const v0 = currentV[idx0];
                        const v1 = currentV[idx1];
                        const v2 = currentV[idx2];
                        
                        // R . v
                        nextV[idx0] = R[0][0]*v0 + R[0][1]*v1 + R[0][2]*v2;
                        nextV[idx1] = R[1][0]*v0 + R[1][1]*v1 + R[1][2]*v2;
                        nextV[idx2] = R[2][0]*v0 + R[2][1]*v1 + R[2][2]*v2;
                    }
                }
                
                // Swap buffers
                const temp = currentV;
                currentV = nextV;
                nextV = temp;
            }
            
            return currentV;
        }
        
        function getTensorLabels(rank) {
            const axes = ['x', 'y', 'z'];
            if (rank === 1) return axes;
            
            let labels = axes;
            for (let r = 2; r <= rank; r++) {
                const newLabels = [];
                for (let l of labels) {
                    for (let a of axes) {
                        newLabels.push(l + a);
                    }
                }
                labels = newLabels;
            }
            return labels;
        }

        // --- Orbit & Ordering Utils ---

        function getLatticeTranslations(groupName) {
            const parts = groupName.split(/\s+/);
            let symbol = "";
            for(let p of parts) {
                if (p.match(/^[A-Z]/) && !p.startsWith("BNS") && !p.startsWith("OG")) {
                    symbol = p;
                    break;
                }
            }
            if (!symbol) {
                const match = groupName.match(/BNS:\s*[\d\.]+\s+([A-Z])/);
                if (match) symbol = match[1];
            }
            
            if (!symbol) return [[0,0,0]];
            
            const firstChar = symbol.charAt(0);
            const translations = [[0,0,0]];
            
            switch(firstChar) {
                case 'I': translations.push([0.5, 0.5, 0.5]); break;
                case 'F': 
                    translations.push([0, 0.5, 0.5]);
                    translations.push([0.5, 0, 0.5]);
                    translations.push([0.5, 0.5, 0]);
                    break;
                case 'A': translations.push([0, 0.5, 0.5]); break;
                case 'B': translations.push([0.5, 0, 0.5]); break;
                case 'C': translations.push([0.5, 0.5, 0]); break;
                case 'R': 
                    translations.push([2/3, 1/3, 1/3]);
                    translations.push([1/3, 2/3, 2/3]);
                    break;
            }
            return translations;
        }

        function getWyckoffOrbit(group, wyckoff, varsIn=null) {
    // Use provided vars (one set for the whole orbit), otherwise make random ones
            const vars = varsIn || { x: Math.random(), y: Math.random(), z: Math.random() };

            
            const parse = (str) => {
                // Simple eval with replacement
                // Note: str might be "1/2-x"
                let s = str.toLowerCase();
                s = s.replace(/x/g, vars.x).replace(/y/g, vars.y).replace(/z/g, vars.z);
                // Handle fractions like 1/2
                // We can use eval() here as data is trusted
                try {
                    return eval(s);
                } catch(e) { return 0; }
            };
            const coordStr = wyckoff.coord;
            const parts = coordStr.split(',');
            const P0 = [parse(parts[0]), parse(parts[1]), parse(parts[2])];
            
            // --- Symbolic Parsing ---
            const parseSymComponent = (s) => {
                // s is like "1/2-x" or "-y"
                // Returns {x:0, y:0, z:0, c:0}
                let res = {x:0, y:0, z:0, c:0};
                s = s.replace(/\s+/g, '').toLowerCase();
                
                // Split by + or - (keeping the sign)
                // Hack: replace - with +- and split by +
                let terms = s.replace(/-/g, '+-').split('+').filter(t => t !== '');
                
                terms.forEach(t => {
                    if (t === '') return;
                    let coeff = 1;
                    if (t.startsWith('-')) { coeff = -1; t = t.substring(1); }
                    if (t.startsWith('+')) { t = t.substring(1); }
                    
                    if (t.includes('x')) { 
                        let val = t.replace('x', '');
                        if (val === '') val = '1';
                        res.x += coeff * parseFloat(val); 
                    }
                    else if (t.includes('y')) { 
                        let val = t.replace('y', '');
                        if (val === '') val = '1';
                        res.y += coeff * parseFloat(val); 
                    }
                    else if (t.includes('z')) { 
                        let val = t.replace('z', '');
                        if (val === '') val = '1';
                        res.z += coeff * parseFloat(val); 
                    }
                    else {
                        // Constant
                        if (t.includes('/')) {
                            const [n, d] = t.split('/');
                            res.c += coeff * (parseFloat(n)/parseFloat(d));
                        } else {
                            res.c += coeff * parseFloat(t);
                        }
                    }
                });
                return res;
            };
            
            const P0_sym = parts.map(parseSymComponent);

            const orbit = [];
            const seen = [];
            
            // Identify Lattice Translations
            const latTrans = getLatticeTranslations(group.name);

            // Apply all group operators AND lattice translations
            latTrans.forEach(lt => {
                group.operators.forEach((op, opIdx) => {
                    // Combine op with lattice translation
                    let tStr = "";
                    if (Math.abs(lt[0])>1e-6 || Math.abs(lt[1])>1e-6 || Math.abs(lt[2])>1e-6) {
                        tStr = ` + (${decimalToFraction(lt[0])},${decimalToFraction(lt[1])},${decimalToFraction(lt[2])})`;
                    }

                    const combinedOp = {
                        R: op.R,
                        t: [op.t[0] + lt[0], op.t[1] + lt[1], op.t[2] + lt[2]],
                        tr: op.tr,
                        str: op.str + tStr
                    };
                    
                    // P' = R*P + t
                    const P_prime = vecAdd(matvec(combinedOp.R, P0), combinedOp.t);
                    
                    // Normalize to [0,1)
                    for(let k=0; k<3; k++) {
                        P_prime[k] = P_prime[k] - Math.floor(P_prime[k] + 1e-6);
                        if (Math.abs(P_prime[k] - 1.0) < 1e-6) P_prime[k] = 0.0;
                    }
                    
                    // Check uniqueness
                    let unique = true;
                    for(let s of seen) {
                        const diff = vecSub(P_prime, s);
                        if (isIntegerVec(diff)) {
                            unique = false;
                            break;
                        }
                    }
                    
                    if (unique) {
                        seen.push(P_prime);
                        
                        // Calculate Symbolic P'
                        // P'_sym = R * P0_sym + t
                        const P_prime_sym = [
                            {x:0, y:0, z:0, c:0},
                            {x:0, y:0, z:0, c:0},
                            {x:0, y:0, z:0, c:0}
                        ];
                        
                        for(let i=0; i<3; i++) {
                            // Row i of R
                            for(let j=0; j<3; j++) {
                                const r = combinedOp.R[i][j];
                                if (r === 0) continue;
                                P_prime_sym[i].x += r * P0_sym[j].x;
                                P_prime_sym[i].y += r * P0_sym[j].y;
                                P_prime_sym[i].z += r * P0_sym[j].z;
                                P_prime_sym[i].c += r * P0_sym[j].c;
                            }
                            P_prime_sym[i].c += combinedOp.t[i];
                            
                            // Normalize constant part to [0,1)
                            let c = P_prime_sym[i].c;
                            c = c - Math.floor(c + 1e-6);
                            if (Math.abs(c - 1.0) < 1e-6) c = 0.0;
                            P_prime_sym[i].c = c;
                        }
                        
                        // Format Symbolic
                    const formatSym = (comp) => {
                        let s = "";
                        // Constant
                        const tol = 1e-4;
                        let c = comp.c;
                        if (Math.abs(c) > tol) {
                            // Fraction
                            const maxDenom = 24;
                            let found = false;
                            for(let d=2; d<=maxDenom; d++) {
                                const num = Math.round(c * d);
                                if (Math.abs(c - num/d) < tol) {
                                    s += `${num}/${d}`;
                                    found = true;
                                    break;
                                }
                            }
                            if (!found) s += c.toFixed(3).replace(/\.?0+$/, '');
                        }
                        
                        // Variables
                        const vars = [
                            {val: comp.x, label: 'x'},
                            {val: comp.y, label: 'y'},
                            {val: comp.z, label: 'z'}
                        ];
                        
                        vars.forEach(v => {
                            if (Math.abs(v.val) > tol) {
                                let sign = (v.val > 0) ? "+" : "-";
                                if (s === "" && sign === "+") sign = "";
                                else if (s !== "" && sign === "+") sign = "+"; // Add + if not first
                                
                                let absVal = Math.abs(v.val);
                                let coeff = "";
                                if (Math.abs(absVal - 1) > tol) {
                                    coeff = absVal.toString(); // e.g. 2x
                                }
                                
                                s += `${sign}${coeff}${v.label}`;
                            }
                        });
                        
                        if (s === "") return "0";
                        return s;
                    };
                    
                    const coordSymStr = P_prime_sym.map(formatSym).join(',');

                    orbit.push({
                        coord: P_prime,
                        coordSym: coordSymStr,
                        op: combinedOp, 
                        opIdx: opIdx
                    });
                }
            });
            });
            
            return orbit;
        }
        function updateOrbitViewer() {
            const wSel = document.getElementById("wyckoffSelect");
            if (!wSel || !currentGroup) return;

            if (wSel.value === "whole") {
                if (orbitGroup3d && scene3d) {
                    scene3d.remove(orbitGroup3d);
                    orbitGroup3d = null;
                }
                clearIsoGroup();
                _lastOrbitForTensor = null;

                // ✅ clear stale sliders too
                clearTensorSlidersUI();
                return;
                }


            const w = currentGroup.wyckoff[parseInt(wSel.value, 10)];
            if (!w) return;

            // If x/y/z are not fixed numbers, pick random values in [0,1)
            // (your wyckoffRandomVars is already meant for this) :contentReference[oaicite:3]{index=3}
            if (!wyckoffRandomVars) {
        wyckoffRandomVars = generateWyckoffRandomVars(currentGroup, w);
            }

            // Build the orbit from the Wyckoff + random vars, then render spheres
            const orbit = getWyckoffOrbit(currentGroup, w, wyckoffRandomVars);
            renderOrbitSpheres(orbit);
            // ✅ Update hover tooltip UP/DOWN/NCL AFTER spheres exist
            const selectedRank = getSelectedRank();
            if (selectedRank) {
            const basis = getBasisForRank(selectedRank);
            if (basis && basis.length > 0) {
                const groups = analyzeAtomOrdering(basis, selectedRank, orbit, currentGroup);
                updateOrbitSphereTooltips(orbit, groups);
            }
            }

            renderTensorIsosurfacesFromSiteTensor(orbit);
            }
        function analyzeAtomOrdering(basis, rank, orbit, group) {
            const dim = Math.pow(3, rank);
            
            // 1. Compute T_i for all atoms for ALL basis vectors
            // To be strictly "Same" or "Opposite", the relationship must hold for the entire general solution space.
            // i.e. for every basis vector b_k, T_i(b_k) == +/- T_j(b_k) with the SAME sign for all k.
            
            const atomTensors = [];
            // Prepare lattice matrices to convert fractional ops to Cartesian
            const M = getIdealLatticeMatrix(group);
            const invM = invertMatrix3x3(M);

            orbit.forEach((atom, i) => {
                const op = atom.op;
                // Use fractional/crystal operator directly for algebra
                const R_frac = op.R;
                const detR = det3x3(R_frac); // determinant is same as Cartesian conjugate
                const tr = (op.tr !== undefined) ? op.tr : 1;
                const factor = (multipoleMode === 'magnetic') ? (tr * detR) : 1;

                const basisTransformed = basis.map(bVec => {
                    let curr = new Float64Array(bVec);
                    let next = new Float64Array(dim);
                    // apply fractional operator directly
                    const transformed = applyOpToTensor(curr, rank, R_frac, next);
                    for(let k=0; k<dim; k++) transformed[k] *= factor;
                    return transformed;
                });
                
                atomTensors.push({
                    id: i,
                    atom: atom,
                    basisVecs: basisTransformed
                });
            });
            
            // 2. Group them
            const groups = [];
            const processed = new Set();
            const TOL = 1e-4;
            
            for(let i=0; i<atomTensors.length; i++) {
                if (processed.has(i)) continue;
                
                const rep = atomTensors[i];
                const group = {
                    id: groups.length,
                    members: []
                };
                
                // Add representative
                group.members.push({
                    atomIdx: i,
                    atom: rep.atom,
                    relation: "Same" // Parallel to group rep
                });
                processed.add(i);
                
                // Find others
                for(let j=i+1; j<atomTensors.length; j++) {
                    if (processed.has(j)) continue;
                    
                    const cand = atomTensors[j];
                    
                    // Check Parallel (Must be parallel for ALL basis vectors)
                    let isParallel = true;
                    for(let b=0; b<basis.length; b++) {
                        const v1 = rep.basisVecs[b];
                        const v2 = cand.basisVecs[b];
                        for(let k=0; k<dim; k++) {
                            if (Math.abs(v1[k] - v2[k]) > TOL) {
                                isParallel = false;
                                break;
                            }
                        }
                        if (!isParallel) break;
                    }
                    
                    if (isParallel) {
                        group.members.push({
                            atomIdx: j,
                            atom: cand.atom,
                            relation: "Same"
                        });
                        processed.add(j);
                        continue;
                    }
                    
                    // Check Antiparallel (Must be antiparallel for ALL basis vectors)
                    let isAntiparallel = true;
                    for(let b=0; b<basis.length; b++) {
                        const v1 = rep.basisVecs[b];
                        const v2 = cand.basisVecs[b];
                        for(let k=0; k<dim; k++) {
                            if (Math.abs(v1[k] + v2[k]) > TOL) {
                                isAntiparallel = false;
                                break;
                            }
                        }
                        if (!isAntiparallel) break;
                    }
                    
                    if (isAntiparallel) {
                        group.members.push({
                            atomIdx: j,
                            atom: cand.atom,
                            relation: "Opposite"
                        });
                        processed.add(j);
                        continue;
                    }
                }
                groups.push(group);
            }
            
            // 3. Also compute relation of every atom to the representative atom (atom 0)
            // This will be used to decide UP/DOWN labels only when tensors are collinear
            const relationToRef = new Array(atomTensors.length).fill("Incomparable");
            if (atomTensors.length > 0) {
                const ref = atomTensors[0];
                for (let j = 0; j < atomTensors.length; j++) {
                    const cand = atomTensors[j];
                    // Compare for Same (equal for all basis vectors) or Opposite (negated for all)
                    let allSame = true;
                    let allOpp = true;
                    for (let b = 0; b < basis.length; b++) {
                        const v1 = ref.basisVecs[b];
                        const v2 = cand.basisVecs[b];
                        for (let k = 0; k < v1.length; k++) {
                            if (Math.abs(v1[k] - v2[k]) > TOL) { allSame = false; break; }
                        }
                        for (let k = 0; k < v1.length; k++) {
                            if (Math.abs(v1[k] + v2[k]) > TOL) { allOpp = false; break; }
                        }
                        if (!allSame && !allOpp) break;
                    }

                    if (allSame) relationToRef[j] = "Same";
                    else if (allOpp) relationToRef[j] = "Opposite";
                    else relationToRef[j] = "Incomparable";
                }
            }

            // If any non-reference atom is Incomparable (non-collinear), then
            // the reference atom should also not be labeled UP/DOWN — mark it Incomparable.
            let anyIncomparable = false;
            for (let j = 1; j < relationToRef.length; j++) {
                if (relationToRef[j] === "Incomparable") { anyIncomparable = true; break; }
            }
            if (anyIncomparable) relationToRef[0] = "Incomparable";

            // Attach relationToRef to group members for convenience
            groups.forEach(g => {
                g.members.forEach(m => {
                    m.relationToRef = relationToRef[m.atomIdx];
                });
            });

            return groups;
        }

        // --- Formatting ---

        function decimalToFraction(val) {
            const tol = 1e-4;
            if (Math.abs(val) < tol) return "0";
            
            const sign = val < -tol ? "-" : "";
            val = Math.abs(val);
            
            // Check for common square roots
            if (Math.abs(val - Math.sqrt(2)) < tol) return sign + "√2";
            if (Math.abs(val - Math.sqrt(3)) < tol) return sign + "√3";
            if (Math.abs(val - 1/Math.sqrt(2)) < tol) return sign + "1/√2";
            if (Math.abs(val - 1/Math.sqrt(3)) < tol) return sign + "1/√3";
            if (Math.abs(val - Math.sqrt(3)/2) < tol) return sign + "√3/2";
            if (Math.abs(val - (Math.sqrt(3) - 1)) < tol) return sign + "(√3-1)";
            if (Math.abs(val - (Math.sqrt(3) + 1)) < tol) return sign + "(√3+1)";
            
            // Check if integer (or close to 1, 2 etc)
            if (Math.abs(val - Math.round(val)) < tol) return sign + Math.round(val).toString();

            const maxDenom = 24; 
            for(let d=2; d<=maxDenom; d++) {
                const num = Math.round(val * d);
                if (Math.abs(val - num/d) < tol) {
                    const gcd = (a, b) => b ? gcd(b, a % b) : a;
                    const common = gcd(num, d);
                    return sign + `${num/common}/${d/common}`;
                }
            }
            return sign + val.toFixed(3).replace(/\.?0+$/, '');
        }

        function displayResult(rank, basis, container, group, wyckoff) {
            const div = document.createElement('div');
            div.className = 'rank-block';
            
            const rankNames = {
                1: "Dipole (Rank 1)",
                2: "Quadrupole (Rank 2)",
                3: "Octupole (Rank 3)",
                4: "Hexadecapole (Rank 4)",
                5: "Triakontadipole (Rank 5)"
            };

            const title = document.createElement('div');
            title.className = 'rank-title';
            title.innerText = rankNames[rank] || `Rank ${rank}`;
            div.appendChild(title);
            // Rank-specific small info (e.g., sublattice translation for Rank 1) will be appended
            // together with the atom ordering block later.
            
            if (basis.length === 0) {
                div.innerHTML += "<div>All components are zero (Forbidden).</div>";
                container.appendChild(div);
                return;
            }

            // Analyze the shape of the general tensor
            const labels = getTensorLabels(rank);
            const nComponents = labels.length;
            const nBasis = basis.length;

            // --- Pre-calculate Component Groups (General Tensor Shape) ---
            // Decide which basis to use for grouping/display: use Cartesian for System results
            let displayBasis = basis;
            if (!wyckoff && group) {
                // Convert each basis vector to Cartesian using the ideal lattice matrix
                try {
                    const latM = getIdealLatticeMatrix(group);
                    const transformTensorToCartesianLocal = (vec, rank, M) => {
                        if (rank === 0) return vec;
                        let out = new Float64Array(vec.length);
                        const n = Math.pow(3, rank);
                        for (let flatOut = 0; flatOut < n; flatOut++) {
                            let idxsOut = [];
                            let rem = flatOut;
                            for (let d = 0; d < rank; d++) {
                                idxsOut.unshift(rem % 3);
                                rem = Math.floor(rem / 3);
                            }
                            let sum = 0;
                            for (let flatIn = 0; flatIn < n; flatIn++) {
                                let idxsIn = [];
                                let remIn = flatIn;
                                for (let d = 0; d < rank; d++) {
                                    idxsIn.unshift(remIn % 3);
                                    remIn = Math.floor(remIn / 3);
                                }
                                let prod = 1;
                                for (let d = 0; d < rank; d++) {
                                    prod *= M[idxsOut[d]][idxsIn[d]];
                                }
                                sum += prod * vec[flatIn];
                            }
                            out[flatOut] = sum;
                        }
                        return out;
                    };

                    displayBasis = basis.map(bv => {
                        try {
                            return transformTensorToCartesianLocal(bv, rank, latM);
                        } catch(e) {
                            console.warn('Cartesian transform failed for system basis', e);
                            return bv;
                        }
                    });
                } catch(e) {
                    console.warn('Failed to compute lattice matrix for system Cartesian display', e);
                    displayBasis = basis;
                }
            }

            // Extract component signatures from the chosen display basis
            const sigs = [];
            for(let i=0; i<nComponents; i++) {
                const s = new Float64Array(nBasis);
                for(let b=0; b<nBasis; b++) s[b] = displayBasis[b][i];
                sigs.push({idx: i, vec: s, norm: Math.sqrt(vecDot(s, s))});
            }
            
            // Group components
            const componentGroups = []; 
            const zeroComponents = [];
            const GROUP_TOL = 1e-4;

            for(let i=0; i<nComponents; i++) {
                const s = sigs[i];
                if (s.norm < GROUP_TOL) {
                    zeroComponents.push(labels[i]);
                    continue;
                }
                
                let found = false;
                for(let g of componentGroups) {
                    const dot = vecDot(s.vec, g.rep.vec);
                    const factor = dot / (g.rep.norm * g.rep.norm);
                    
                    let resSq = 0;
                    for(let k=0; k<nBasis; k++) {
                        const diff = s.vec[k] - factor * g.rep.vec[k];
                        resSq += diff * diff;
                    }
                    
                    if (Math.sqrt(resSq) < GROUP_TOL) {
                        g.members.push({idx: i, factor: factor});
                        found = true;
                        break;
                    }
                }
                
                if (!found) {
                    componentGroups.push({
                        rep: s,
                        members: [{idx: i, factor: 1.0}]
                    });
                }
            }
            
            // Sort members within groups
            componentGroups.forEach(g => {
                g.members.sort((a, b) => {
                    const aOne = Math.abs(a.factor - 1.0) < GROUP_TOL;
                    const bOne = Math.abs(b.factor - 1.0) < GROUP_TOL;
                    if (aOne && !bOne) return -1;
                    if (!aOne && bOne) return 1;
                    const aMinus = Math.abs(a.factor + 1.0) < GROUP_TOL;
                    const bMinus = Math.abs(b.factor + 1.0) < GROUP_TOL;
                    if (aMinus && !bMinus) return -1;
                    if (!aMinus && bMinus) return 1;
                    return a.idx - b.idx;
                });
            });

            // --- Atom Ordering Calculation ---
            if (group && wyckoff) {
                const orderingDiv = document.createElement('div');
                orderingDiv.style.marginBottom = "15px";
                orderingDiv.style.padding = "10px";
                orderingDiv.style.backgroundColor = "#f9f9f9";
                orderingDiv.style.border = "1px solid #eee";
                
                orderingDiv.innerHTML = "<strong>Atom Ordering (Multipole Relations):</strong><br>";
                
                // 1. Generate Orbit
                const orbit = _lastOrbitForTensor || getWyckoffOrbit(group, wyckoff, wyckoffRandomVars);

                
                    // 2. Calculate relations
                if (basis.length > 0) {
                    const groups = analyzeAtomOrdering(basis, rank, orbit, group);
                    updateOrbitSphereTooltips(orbit, groups);

                        // If rank==1, compute sublattice translation between representative and an Opposite member
                        if (rank === 1) {
                            let sublatticeInfo = null;
                            // Determine if any incomparable present
                            let anyIncomparableLocal = false;
                            for (let g of groups) {
                                g.members.forEach(m => { if (m.relationToRef === 'Incomparable') anyIncomparableLocal = true; });
                            }

                            if (!anyIncomparableLocal) {
                                // Build lists of up (Same) and down (Opposite) atom indices
                                const upIdx = [];
                                const downIdx = [];
                                for (let g of groups) {
                                    g.members.forEach(m => {
                                        if (m.relationToRef === 'Same') upIdx.push(m.atomIdx);
                                        else if (m.relationToRef === 'Opposite') downIdx.push(m.atomIdx);
                                    });
                                }

                                // If we have both up and down atoms, search for a pure translation t
                                // (R == identity) among group operators combined with lattice translations
                                if (upIdx.length > 0 && downIdx.length > 0) {
                                    const latTrans = getLatticeTranslations(group.name);
                                    const identityR = [[1,0,0],[0,1,0],[0,0,1]];
                                    const isIdentity = (R) => {
                                        const tol = 1e-6;
                                        for (let i=0;i<3;i++) for (let j=0;j<3;j++) if (Math.abs(R[i][j] - identityR[i][j]) > tol) return false;
                                        return true;
                                    };

                                    const coords = orbit.map(o => o.coord);

                                    // helper to normalize vector to [0,1)
                                    const normFrac = (v) => v.map((c)=>{ let x = c - Math.floor(c + 1e-6); if (Math.abs(x-1.0) < 1e-6) x = 0.0; return x; });
                                    const approxEqual = (a,b) => { const tol=1e-4; for(let k=0;k<3;k++) if (Math.abs(a[k]-b[k])>tol) return false; return true; };

                                    outer: for (let lt of latTrans) {
                                        for (let op of group.operators) {
                                            // combined op
                                            const t = [op.t[0] + lt[0], op.t[1] + lt[1], op.t[2] + lt[2]];
                                            if (!isIdentity(op.R)) continue; // require pure translation

                                            const tnorm = normFrac(t);

                                            // Check that applying t to each up atom lands on some down atom
                                            let allMap = true;
                                            for (let ui of upIdx) {
                                                const target = normFrac([coords[ui][0] + tnorm[0], coords[ui][1] + tnorm[1], coords[ui][2] + tnorm[2]]);
                                                let found = false;
                                                for (let di of downIdx) {
                                                    if (approxEqual(target, coords[di])) { found = true; break; }
                                                }
                                                if (!found) { allMap = false; break; }
                                            }

                                            if (allMap) {
                                                sublatticeInfo = tnorm;
                                                break outer;
                                            }
                                        }
                                    }
                                }
                            }

                            // Attach sublatticeInfo to orderingDiv for display later
                            orderingDiv.sublatticeInfo = sublatticeInfo; // may be null
                            orderingDiv.anyIncomparableLocal = anyIncomparableLocal;
                        }
                    
                    const grid = document.createElement('div');
                    grid.style.display = "grid";
                    grid.style.gridTemplateColumns = "repeat(auto-fill, minmax(280px, 1fr))";
                    grid.style.gap = "10px";
                    grid.style.marginTop = "10px";
                    
                    // Flatten for display order (by atom index)
                    const flatList = new Array(orbit.length);
                    groups.forEach(g => {
                        // Determine group color: green when all Same, red when any Opposite (and no Incomparable), neutral otherwise
                        let gColor = 'neutral';
                        const rels = g.members.map(mm => mm.relation || mm.relationToRef);
                        if (rels.some(r => r === 'Incomparable')) gColor = 'neutral';
                        else if (rels.every(r => r === 'Same')) gColor = 'green';
                        else if (rels.some(r => r === 'Opposite')) gColor = 'red';

                        g.members.forEach(m => {
                            let bgCol, borderCol, textCol;
                            
                            // Decide display color based on relation to reference atom (atom 1)
                            // Only show UP/DOWN when tensors are collinear (Same or Opposite). Otherwise show neutral styling.
                            const relToRef = m.relationToRef || m.relation; // fallback
                            if (relToRef === "Same") {
                                // Red for UP
                                bgCol = "#ffebee"; // Very light red
                                borderCol = "#ef5350"; // Red
                                textCol = "#b71c1c"; // Dark red
                            } else if (relToRef === "Opposite") {
                                // Blue for DOWN
                                bgCol = "#e3f2fd"; // Very light blue
                                borderCol = "#42a5f5"; // Blue
                                textCol = "#0d47a1"; // Dark blue
                            } else {
                                // Neutral styling when not collinear
                                bgCol = "#f5f5f5";
                                borderCol = "#cccccc";
                                textCol = "#555555";
                            }

                            flatList[m.atomIdx] = {
                                member: m,
                                groupId: g.id,
                                groupSize: g.members.length,
                                groupColor: gColor,
                                style: `background-color: ${bgCol}; color: ${textCol}; border: 1px solid ${borderCol};`
                            };
                        });
                    });

                    // If any atom in the orbit is non-collinear relative to the reference,
                    // color all atoms with the non-collinear (neutral) styling per request.
                    let anyIncomparable = false;
                    for (let i = 0; i < flatList.length; i++) {
                        const it = flatList[i];
                        if (!it) continue;
                        const rel = it.member.relationToRef || it.member.relation;
                        if (rel === "Incomparable") { anyIncomparable = true; break; }
                    }
                    if (anyIncomparable) {
                        for (let i = 0; i < flatList.length; i++) {
                            if (!flatList[i]) continue;
                            flatList[i].style = `background-color: #f5f5f5; color: #555555; border: 1px solid #cccccc;`;
                        }
                    }
                    
                    // Container for selected tensor
          // Container for selected tensor will be created on demand when an atom is clicked
            let selectedDiv = null;
            // Variable registry for this ordering run: maps canonical variable label (e.g. 'x')
            // to an object { label: "x'", origCoefNum: 1.0, origCoefStr: "1" , defined: true }
            // This lets later atoms reuse the same short label and express coefficients
            // relative to the original coefficient rather than creating new primed labels.
            // use global VARIABLE_REGISTRY instead of a local one

          const ensureSelectedDiv = () => {
            if (!selectedDiv) {
              selectedDiv = document.createElement('div');
              selectedDiv.id = `selected-tensor-rank-${rank}`;
              selectedDiv.style.marginTop = "15px";
              selectedDiv.style.padding = "10px";
              selectedDiv.style.border = "1px solid #ccc";
              selectedDiv.style.backgroundColor = "#fff";
              selectedDiv.style.display = "none"; // Hidden initially
              orderingDiv.appendChild(selectedDiv);
            }
          };
                    
                    const showTensorForAtom = (atom) => {
            try {
            const dim = Math.pow(3, rank);
                        const op = atom.op;
                        // Prepare lattice matrices
                        const latM = getIdealLatticeMatrix(group);
                        const invLatM = invertMatrix3x3(latM);
                        // Work in fractional (reduced) coordinates for the algebra
                        const R_frac = op.R;
                        const detR = det3x3(R_frac);
                        const tr = (op.tr !== undefined) ? op.tr : 1;
                        const factor = (multipoleMode === 'magnetic') ? (tr * detR) : 1;

                        // 1. Calculate transformed basis vectors in fractional (reduced) coordinates
                        let transformedBasisFrac = basis.map(bVec => {
                            let curr = new Float64Array(bVec);
                            let next = new Float64Array(dim);
                            const transformed = applyOpToTensor(curr, rank, R_frac, next);
                            for(let k=0; k<dim; k++) transformed[k] *= factor;
                            return transformed;
                        });

                        

                        // 2. Optionally convert to display basis
                        const showCrystal = document.getElementById('crystalBasisSwitch')?.checked;
                        // helper: transform tensor components by a 3x3 matrix M (applies M on each index)
                        function transformTensorToCartesian(vec, rank, M) {
                            if (rank === 0) return vec;
                            let out = new Float64Array(vec.length);
                            const n = Math.pow(3, rank);
                            for (let flatOut = 0; flatOut < n; flatOut++) {
                                let idxsOut = [];
                                let rem = flatOut;
                                for (let d = 0; d < rank; d++) {
                                    idxsOut.unshift(rem % 3);
                                    rem = Math.floor(rem / 3);
                                }
                                let sum = 0;
                                for (let flatIn = 0; flatIn < n; flatIn++) {
                                    let idxsIn = [];
                                    let remIn = flatIn;
                                    for (let d = 0; d < rank; d++) {
                                        idxsIn.unshift(remIn % 3);
                                        remIn = Math.floor(remIn / 3);
                                    }
                                    let prod = 1;
                                    for (let d = 0; d < rank; d++) {
                                        prod *= M[idxsOut[d]][idxsIn[d]];
                                    }
                                    sum += prod * vec[flatIn];
                                }
                                out[flatOut] = sum;
                            }
                            return out;
                        }

                        // Build Cartesian-display versions by converting the fractional basis
                        // (we keep the algebra in reduced coordinates and convert only for display)
                        const transformedBasisCart = transformedBasisFrac.map(bVec => transformTensorToCartesian(bVec, rank, latM));
                        // Decide which basis to display: fractional (crystal) or Cartesian
                        const transformedBasis = showCrystal ? transformedBasisFrac : transformedBasisCart;

                        // If every transformed basis vector is (near) zero, mark as forbidden
                        const ZERO_TOL = 1e-6; // relax tolerance for higher-rank numeric noise
                        let allZero = true;
                        // Check both Cartesian-transformed basis and display basis to be robust
                        for (let b = 0; b < transformedBasis.length; b++) {
                            const vecDisp = transformedBasis[b];
                            const vecCart = transformedBasisCart[b] || transformedBasis[b];
                            let normDisp = 0, normCart = 0;
                            for (let i = 0; i < vecDisp.length; i++) {
                                normDisp += vecDisp[i] * vecDisp[i];
                                normCart += vecCart[i] * vecCart[i];
                            }
                            if (Math.sqrt(normDisp) > ZERO_TOL || Math.sqrt(normCart) > ZERO_TOL) { allZero = false; break; }
                        }
                        if (allZero) {
              ensureSelectedDiv();
              selectedDiv.innerHTML = '<div>All components are zero (Forbidden).</div>';
              selectedDiv.style.display = 'block';
                            return;
                        }
        // Crystal/Cartesian basis switch event
        const crystalBasisSwitch = document.getElementById('crystalBasisSwitch');
        if (crystalBasisSwitch) {
            crystalBasisSwitch.onchange = triggerCompute;
        }

                        // 2. Map basis coefficients to component labels
                        // We have N independent variables (componentGroups.length)
                        // Variable k corresponds to componentGroups[k].members[0].idx (the representative)
                        // Let's call the variable name L_k (e.g. "x", "Qxy")
                        // We need to express the atom tensor components in terms of L_k.
                        
                        // Construct Matrix M where M_kb = basis[b][rep_idx_k]
                        // This maps basis coeffs c to variable values v: v = M * c
                        const N = componentGroups.length;
                        const M = [];
                        for(let k=0; k<N; k++) {
                            const repIdx = componentGroups[k].members[0].idx;
                            const row = [];
                            for(let b=0; b<nBasis; b++) {
                                row.push(basis[b][repIdx]);
                            }
                            M.push(row);
                        }
                        
            // Invert/pseudo-inverse M to get c = M_inv * v
            // M has shape N x nBasis (N variables rows, nBasis columns). It may not be square.
            // Compute Moore-Penrose style pseudoinverse using closed-form depending on ranks:
            // - If N <= nBasis: pinv = M^T * inv(M * M^T)
            // - If N > nBasis:  pinv = inv(M^T * M) * M^T
            const transpose = (A) => A[0].map((_, c) => A.map(r => r[c]));
            const matMul = (A, B) => {
              const p = A.length, q = A[0].length, r = B[0].length;
              const out = Array.from({length: p}, () => Array(r).fill(0));
              for (let i=0;i<p;i++) for (let k=0;k<q;k++) if (A[i][k] !== 0) for (let j=0;j<r;j++) out[i][j] += A[i][k] * B[k][j];
              return out;
            };

            const addDiag = (A, eps) => A.map((row, i) => row.map((v, j) => v + (i===j?eps:0)));

            const invertSquare = (A) => {
              // Reuse local Gaussian-elimination based inverter (works for reasonably sized square matrices)
              const n = A.length;
              const aug = A.map((row, i) => [...row, ...Array(n).fill(0).map((_, j) => i===j?1:0)]);
              for(let i=0;i<n;i++){
                let pivot = aug[i][i];
                if (Math.abs(pivot) < 1e-12) {
                  let swapped = false;
                  for(let j=i+1;j<n;j++){
                    if (Math.abs(aug[j][i]) > 1e-12) { [aug[i], aug[j]] = [aug[j], aug[i]]; pivot = aug[i][i]; swapped = true; break; }
                  }
                  if (!swapped) throw new Error('Singular matrix');
                }
                for(let j=0;j<2*n;j++) aug[i][j] /= pivot;
                for(let k=0;k<n;k++) if(k!==i){ const f = aug[k][i]; for(let j=0;j<2*n;j++) aug[k][j] -= f * aug[i][j]; }
              }
              return aug.map(row => row.slice(n));
            };

            const pseudoInverse = (Mmat) => {
              const rows = Mmat.length; const cols = Mmat[0].length;
              const MT = transpose(Mmat);
              // If square try direct inversion first
              if (rows === cols) {
                try { return invertSquare(Mmat); } catch(e) { /* fall back */ }
              }
              if (rows <= cols) {
                // pinv = MT * inv(M * MT)
                let A = matMul(Mmat, MT); // rows x rows
                try {
                  return matMul(MT, invertSquare(A));
                } catch(e) {
                  // regularize
                  const eps = 1e-8;
                  A = addDiag(A, eps);
                  return matMul(MT, invertSquare(A));
                }
              } else {
                // rows > cols -> pinv = inv(MT * M) * MT
                let A = matMul(MT, Mmat); // cols x cols
                try {
                  return matMul(invertSquare(A), MT);
                } catch(e) {
                  const eps = 1e-8;
                  A = addDiag(A, eps);
                  return matMul(invertSquare(A), MT);
                }
              }
            };

            let M_inv;
            try {
              M_inv = pseudoInverse(M);
            } catch(e) {
              console.error('Failed to compute pseudoinverse of M', e);
              return;
            }
                        
                        // Calculate Mode Vectors W_k for each variable L_k
                        // W_k = sum_b (M_inv[b][k] * TransformedBasis[b])
                        const W = [];
                        for(let k=0; k<N; k++) {
                            const vec = new Float64Array(dim);
                            for(let b=0; b<nBasis; b++) {
                                const coeff = M_inv[b][k];
                                for(let d=0; d<dim; d++) {
                                    vec[d] += coeff * transformedBasis[b][d];
                                }
                            }
                            W.push(vec);
                        }

                        // Also compute Cartesian-mode vectors for 3D display using the same coefficients
                        const W_cart = [];
                        for(let k=0; k<N; k++) {
                            const vec = new Float64Array(dim);
                            for(let b=0; b<nBasis; b++) {
                                const coeff = M_inv[b][k];
                                for(let d=0; d<dim; d++) {
                                    vec[d] += coeff * transformedBasisCart[b][d];
                                }
                            }
                            W_cart.push(vec);
                        }
                        
                        // Now construct the display strings
                        const comps = new Array(dim).fill("0");
                        const colors = new Array(dim).fill("#999");
                        
                        // Generate colors for variables (same as main display)
                        const varColors = componentGroups.map((_, i) => `hsl(${(i * 137.508) % 360}, 80%, 35%)`);
                        
            // Build parts for each component and keep them for later use
            const origParts = new Array(dim).fill(null);
            for(let d=0; d<dim; d++) {
              let parts = [];
              for(let k=0; k<N; k++) {
                const val = W[k][d];
                if (Math.abs(val) > 1e-5) {
                  const repIdx = componentGroups[k].members[0].idx;
                  const label = labels[repIdx];
                  parts.push({val: val, label: label, colIdx: k});
                }
              }
              origParts[d] = parts;
                            
              if (parts.length > 0) {
                let s = "";
                parts.forEach((p, idx) => {
                  let val = p.val;
                  let sign = (val >= 0) ? "+" : "-";
                                    
                  if (idx === 0) {
                    if (sign === "-") s += "-";
                  } else {
                    s += (sign === "+") ? " + " : " - ";
                  }
                                    
                  let absVal = Math.abs(val);
                  let valStr = "";
                  if (Math.abs(absVal - 1) < 1e-5) {
                    valStr = "";
                  } else {
                    let frac = decimalToFraction(absVal);
                    if (frac.startsWith("-")) frac = frac.substring(1);
                    valStr = "(" + frac + ")";
                  }
                                    
                  s += valStr + p.label;
                });
                comps[d] = s;
                                
                if (parts.length === 1) {
                  colors[d] = varColors[parts[0].colIdx];
                } else {
                  colors[d] = "#000";
                }
              }
            }
                        
            // Replace repeated expressions with short labels derived from component indices.
            // Improved behavior: if a component is a single variable times a scalar
            // (e.g. "-1/2 x"), group it together with the bare variable ("x") so that
            // the label def becomes x' = x and other components are written as e.g. -1/2 x'.
            

        const axes = ['x','y','z'];
        const exprGroups = {}; // key -> list of component indices

        // Build variable label -> variable index map (k -> labels)
        const varLabelToK = {};
        for (let k = 0; k < N; k++) {
          const repIdx = componentGroups[k].members[0].idx;
          varLabelToK[ labels[repIdx] ] = k;
        }

        // Detect proportional variable modes using W_cart (Cartesian mode vectors)
        const varVecs = W_cart; // array of Float64Array length N
        const VAR_EQ_TOL = 1e-6;
        const parent = new Array(N);
        for (let i=0;i<N;i++) parent[i]=i;
        function find(x){ while(parent[x]!==x){ parent[x]=parent[parent[x]]; x=parent[x]; } return x; }
        function union(a,b){ const pa=find(a), pb=find(b); if(pa!==pb) parent[pb]=pa; }

        const dot = (a,b) => { let s=0; for (let i=0;i<a.length;i++) s+=a[i]*b[i]; return s; };
        const norm = (a) => Math.sqrt(dot(a,a));

        for (let i=0;i<N;i++) {
          const vi = varVecs[i];
          const ni = norm(vi);
          if (ni < 1e-9) continue;
          for (let j=i+1;j<N;j++) {
            const vj = varVecs[j];
            const nj = norm(vj);
            if (nj < 1e-9) continue;
            const c = dot(vj, vi) / dot(vi, vi);
            // residual || vj - c*vi ||
            let res = 0;
            for (let m=0;m<vi.length;m++) {
              const diff = vj[m] - c*vi[m];
              res += diff*diff;
            }
            res = Math.sqrt(res);
            if (res < VAR_EQ_TOL * ni) {
              union(i,j);
            }
          }
        }

        // compute canonical representative and scale to rep
        const repMap = new Array(N).fill(0);
        for (let k=0;k<N;k++) repMap[k] = find(k);
        const groupReps = {};
        for (let k=0;k<N;k++) {
          const r = repMap[k];
          if (!groupReps[r]) groupReps[r] = [];
          groupReps[r].push(k);
        }
        const scaleToRep = new Array(N).fill(1.0);
        for (const rStr of Object.keys(groupReps)) {
          const members = groupReps[rStr];
          // choose smallest index as canonical rep
          members.sort((a,b)=>a-b);
          const rep = members[0];
          const vre = varVecs[rep];
          const denom = dot(vre, vre);
          for (const m of members) {
            if (m === rep) { scaleToRep[m] = 1.0; continue; }
            const c = denom > 1e-12 ? dot(varVecs[m], vre) / denom : 1.0;
            scaleToRep[m] = c;
          }
        }

        // Now build exprGroups using canonical labels; for single-variable parts, group by
        // the canonical variable only (not the numeric magnitude). The numeric magnitude
        // is stored in exprMeta as the original coefficient seen for that canonical
        // variable (the first occurrence) and will be used when formatting entries
        // so later atoms can reuse the same short label as a Global Reference Variable.
        // We keep an exprMeta lookup for repLabel keys and a compCoeff array to
        // preserve per-component numeric coefficients.
        const exprMeta = {};
        const compCoeff = new Array(comps.length).fill(0);
        for (let d = 0; d < comps.length; d++) {
          const e = comps[d];
          if (e === '0') continue;
          const parts = origParts[d] || [];
          if (parts.length === 1) {
            const k = parts[0].colIdx;
            const rep = repMap[k];
            const repLabel = labels[ componentGroups[rep].members[0].idx ];
            const coeff = parts[0].val * (scaleToRep[k] || 1.0);
            const absFrac = decimalToFraction(Math.abs(coeff));
            const key = repLabel; // group by the canonical variable only
            if (!exprGroups[key]) exprGroups[key] = [];
            exprGroups[key].push(d);
            // Record the FIRST observed coefficient for this canonical variable as
            // the canonical definition coefficient (origCoef). Do not overwrite on
            // subsequent observations so the first label stays the global one.
            if (!exprMeta[key]) exprMeta[key] = { repLabel: repLabel, origCoefStr: absFrac, origCoefNum: Math.abs(coeff) };
            compCoeff[d] = coeff;
          } else {
            if (!exprGroups[e]) exprGroups[e] = [];
            exprGroups[e].push(d);
          }
        }
            const usedLabels = {};
            const compDefs = [];
            const labelMap = {};

            const idxsFromFlat = (flat, r) => {
              const idxs = [];
              let tmp = flat;
              for (let n = 0; n < r; n++) { idxs.unshift(tmp % 3); tmp = Math.floor(tmp / 3); }
              return idxs;
            };

            for (const [key, list] of Object.entries(exprGroups)) {
              // If the key corresponds to a canonical variable (single-var groups),
              // prefer to use the component index of that canonical variable as the
              // base for the generated short label. Otherwise fall back to the
              // smallest component index found in the group.
              let repCompIdx;
              if (varLabelToK[key] !== undefined) {
                const varK = varLabelToK[key];
                // component index for this variable's representative
                repCompIdx = componentGroups[varK].members[0].idx;
              } else {
                repCompIdx = Math.min(...list);
              }

              const idxs = idxsFromFlat(repCompIdx, rank);
              let base = idxs.map(ii => axes[ii]).join('');
              if (base.length < rank) base = base + base[base.length-1].repeat(rank - base.length);
              // If this key is a canonical variable, prefer the bare base label
              // (no prime) so coefficients display as e.g. (√3/2)zxxyz rather than
              // (√3/2)zxxyz'. Otherwise use the primed short label. However,
              // if a GLOBAL variableRegistry already holds a label for this
              // canonical variable, reuse it instead of inventing a new one.
              let label = (varLabelToK[key] !== undefined) ? base : base + "'";

              // If the variableRegistry already has a label for this canonical
              // variable, prefer that label and avoid creating a new primed name.
              if (varLabelToK[key] !== undefined && VARIABLE_REGISTRY[key] && VARIABLE_REGISTRY[key].label) {
                label = VARIABLE_REGISTRY[key].label;
              } else {
                // Ensure local uniqueness (within this atom) and avoid collisions
                // with labels already registered in variableRegistry for OTHER
                // canonical variables.
                const labelInUseGlobally = Object.keys(VARIABLE_REGISTRY).find(r => VARIABLE_REGISTRY[r].label === label && r !== key);
                if (labelInUseGlobally) {
                  let k = 1;
                  let newLabel = `${label}_${k}`;
                  while ((usedLabels[newLabel] && usedLabels[newLabel] !== key) || Object.keys(VARIABLE_REGISTRY).some(r => VARIABLE_REGISTRY[r].label === newLabel)) { k++; newLabel = `${label}_${k}`; }
                  label = newLabel;
                } else if (usedLabels[label] && usedLabels[label] !== key) {
                  let k = 1;
                  let newLabel = `${label}_${k}`;
                  while (usedLabels[newLabel] && usedLabels[newLabel] !== key) { k++; newLabel = `${label}_${k}`; }
                  label = newLabel;
                }
              }

                // Decide the definition expression for the label. If this exprGroup
                // corresponds to a single-variable with numeric magnitude (we stored
                // it in exprMeta), move the numeric part into the definition and
                // keep the tensor entries coefficient-free (sign stays with entry).
                const meta = exprMeta[key];
                if (meta) {
                  // meta.repLabel is the canonical variable, meta.origCoefStr is the
                  // magnitude observed when this label was first created. If this
                  // canonical variable is already registered globally, do NOT
                  // re-add the definition locally (avoid duplicate defs). Otherwise
                  // record the definition (moving magnitude to the definition).
                  if (VARIABLE_REGISTRY[key] && VARIABLE_REGISTRY[key].label === label && VARIABLE_REGISTRY[key].defined) {
                    // already globally defined; nothing to add locally
                  } else {
                    if (meta.origCoefStr === '1') usedLabels[label] = meta.repLabel;
                    else usedLabels[label] = `(${meta.origCoefStr})` + meta.repLabel;
                    // Register globally so later atoms will reuse this label
                    VARIABLE_REGISTRY[key] = VARIABLE_REGISTRY[key] || {};
                    VARIABLE_REGISTRY[key].label = label;
                    VARIABLE_REGISTRY[key].origCoefNum = meta.origCoefNum || 1.0;
                    VARIABLE_REGISTRY[key].origCoefStr = meta.origCoefStr || '1';
                    VARIABLE_REGISTRY[key].defined = true;
                  }
                } else {
                  if (varLabelToK[key] !== undefined) {
                      // For canonical variables with no simple numeric magnitude,
                      // register/reuse a label if one exists; otherwise create a
                      // local definition mapping the generated short label to the
                      // canonical variable name.
                      if (VARIABLE_REGISTRY[key] && VARIABLE_REGISTRY[key].label) {
                        usedLabels[ VARIABLE_REGISTRY[key].label ] = key;
                        label = VARIABLE_REGISTRY[key].label;
                      } else {
                        usedLabels[label] = key;
                        VARIABLE_REGISTRY[key] = VARIABLE_REGISTRY[key] || {};
                        VARIABLE_REGISTRY[key].label = label;
                        VARIABLE_REGISTRY[key].origCoefNum = 1.0;
                        VARIABLE_REGISTRY[key].origCoefStr = '1';
                        VARIABLE_REGISTRY[key].defined = true;
                      }
                  } else {
                    // Use the full expression string from the representative component
                    usedLabels[label] = comps[repCompIdx];
                  }
                }

                  // Map canonical repLabel -> chosen short label (used for lookup
                  // when rendering entries).  Note: we key by repLabel (not by
                  // coefficient) so that different numeric multiples reuse the same
                  // short label.
                  labelMap[key] = label;
            }

            // Replace components: for single-var keys remove numeric magnitude from
            // the tensor entry (keep sign only) and move the magnitude into the
            // label definition (created above). Multi-term groups still map to
            // their short label if available.
            for (let d = 0; d < comps.length; d++) {
              const e = comps[d];
              if (e === '0') continue;
              const parts = origParts[d] || [];
              if (parts.length === 1) {
                const k = parts[0].colIdx;
                const rep = repMap[k];
                const repLabel = labels[ componentGroups[rep].members[0].idx ];
                const coeff = compCoeff[d] || parts[0].val * (scaleToRep[k] || 1.0);
                const absVal = Math.abs(coeff);
                const coefStr = decimalToFraction(absVal);
                // Lookup the short label by the canonical repLabel (ignore coefficient)
                const lbl = labelMap[repLabel] || (VARIABLE_REGISTRY[repLabel] && VARIABLE_REGISTRY[repLabel].label);
                if (!lbl) continue;

                // Compute the displayed multiplicative factor relative to the
                // originally registered coefficient for this canonical variable.
                const origCoef = (VARIABLE_REGISTRY[repLabel] && VARIABLE_REGISTRY[repLabel].origCoefNum) || (exprMeta[repLabel] && exprMeta[repLabel].origCoefNum) || 1.0;
                const ratio = coeff / (origCoef || 1.0);
                const sign = ratio < 0 ? '-' : '';
                const absRatio = Math.abs(ratio);
                let ratioStr = '';
                if (Math.abs(absRatio - 1) > 1e-8) {
                  ratioStr = decimalToFraction(absRatio);
                  if (ratioStr === '1') ratioStr = '';
                  else ratioStr = `${ratioStr}`;
                }

                // If ratioStr is empty (i.e., 1), show just the label (with sign).
                comps[d] = `${sign}${ratioStr}${lbl}`;
                // Use canonical variable color
                colors[d] = varColors[rep];
              } else {
                // multi-term: replace full expression with label if available
                const lbl = labelMap[e];
                if (lbl) comps[d] = lbl;
              }
            }

            for (let lbl in usedLabels) {
              // Omit trivial self-definitions like `zxxyz = zxxyz` to avoid clutter
              if (usedLabels[lbl] === lbl) continue;
              compDefs.push({lbl: lbl, expr: usedLabels[lbl]});
            }

            // Assign consistent colors to entries that share the same replacement label
            const labelColorMap = {};
            for (let i = 0; i < compDefs.length; i++) {
              const hue = (i * 137.508) % 360; // generate distinct hues
              labelColorMap[compDefs[i].lbl] = `hsl(${hue}, 80%, 35%)`;
            }
            for (let d = 0; d < comps.length; d++) {
              const e = comps[d];
              if (e !== '0' && labelColorMap[e]) colors[d] = labelColorMap[e];
            }

            // Format
                        const coordStr = atom.coordSym || atom.coord.map(x => decimalToFraction(x)).join(',');
                        let html = `<h4>Tensor Shape for Atom at ${coordStr}</h4>`;
                        html += `<div style="font-size: 0.9em; margin-bottom: 5px;">
                                    <strong>Variables:</strong> <span style="font-family: monospace;">
                                    ${componentGroups.map((g, i) => `<span style="color:${varColors[i]}">${labels[g.members[0].idx]}</span>`).join(', ')}
                                    </span>
                                 </div>`;
                        html += `<pre style="font-family: monospace; background: #eee; padding: 10px; overflow-x: auto;">`;
                        
                        const useVoigt = document.getElementById('voigtCheck').checked;
                        
                        if (useVoigt) {
                            if (rank === 1) html += formatSymbolicMatrix1(comps, colors);
                            else if (rank === 2) html += formatVoigtMatrix2(comps, colors);
                            else if (rank === 3) html += formatVoigtMatrix3(comps, colors);
                            else if (rank === 4) html += formatVoigtMatrix4(comps, colors);
                            else if (rank === 5) html += formatVoigtMatrix5(comps, colors);
                        } else {
                            if (rank === 1) html += formatSymbolicMatrix1(comps, colors);
                            else if (rank === 2) html += formatSymbolicMatrix2(comps, colors);
                            else if (rank === 3) html += formatSymbolicMatrix3(comps, colors);
                            else if (rank === 4) html += formatSymbolicMatrix4(comps, colors);
                            else if (rank === 5) html += formatSymbolicMatrix5(comps, colors);
                        }
                        
                        html += `</pre>`;
            if (compDefs.length > 0) {
              // Add label definitions block (hidden by default). Visibility is controlled by header toggle (#showLabelDefsSwitch)
              html += `<div class="label-defs" style="display:none; font-size:0.9em; margin-top:6px;"><strong>Label definitions:</strong><br><pre style="font-family: monospace; background: #f6f6f6; padding:8px;">`;
              compDefs.forEach(d => { html += `${d.lbl} = ${d.expr}\n`; });
              html += `</pre></div>`;
            }
                        
      ensureSelectedDiv();
      selectedDiv.innerHTML = html;
      // Respect the header toggle state immediately when a new selectedDiv is created
      const headerToggle = document.getElementById('showLabelDefsSwitch');
      if (headerToggle && headerToggle.checked) {
        const defs = selectedDiv.querySelectorAll('.label-defs');
        defs.forEach(d => d.style.display = 'block');
      }
      selectedDiv.style.display = "block";
            } catch (err) {
              console.error('Error in showTensorForAtom:', err);
              ensureSelectedDiv();
              selectedDiv.innerHTML = `<div style="color:red;font-weight:700;">Error rendering tensor: ${err && err.message ? err.message : err}</div>`;
              selectedDiv.style.display = 'block';
            }
          };

                    for(let i=0; i<orbit.length; i++) {
                        const item = flatList[i];
                        if (!item) continue; 
                        
                        const m = item.member;
                        const p = m.atom;
                        const coordStr = p.coordSym || p.coord.map(x => decimalToFraction(x)).join(',');
                        const opStr = p.op.str || "Unknown";
                        
                        let relText = "";
                        const relToRef = m.relationToRef || m.relation;
                        if (anyIncomparable) {
                            relText = "NCL"; // Force neutral (Non-ColLinear) when any non-collinear present
                        } else {
                            if (relToRef === "Same") relText = "UP (+)";
                            else if (relToRef === "Opposite") relText = "DOWN (-)";
                            else relText = "NCL"; // Not collinear / no UP/DOWN
                        }
                        
                        let atomLabel = `Atom ${i+1}`;
                        if (i === 0) atomLabel += " (Rep.)";

                        const card = document.createElement('div');
                        card.style.cssText = `padding: 8px; border-radius: 4px; cursor: pointer; transition: transform 0.1s; ${item.style}`;
                        card.onmouseover = () => card.style.transform = "scale(1.02)";
                        card.onmouseout = () => card.style.transform = "scale(1)";
                        card.onclick = () => showTensorForAtom(p);

                        card.innerHTML = `
                                <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 4px;">
                                    <span style="font-weight: bold;">${atomLabel}</span>
                                    ${ item.groupSize && item.groupSize > 1 ? (() => {
                                        let badgeStyle = `font-size: 0.85em; font-weight: bold; padding: 2px 6px; border-radius: 10px; border: 1px solid rgba(0,0,0,0.1);`;
                                        if (item.groupColor === 'green') badgeStyle += ' background: #e8f5e9; color: #1b5e20; border-color: #a5d6a7;';
                                        else if (item.groupColor === 'red') badgeStyle += ' background: #ffebee; color: #b71c1c; border-color: #ef9a9a;';
                                        else badgeStyle += ' background: rgba(255,255,255,0.6); color: #000;';
                                        return `<span style="${badgeStyle}">Group ${item.groupId + 1}</span>`;
                                    })() : '' }
                                </div>
                                <div style="font-size: 0.9em;">Pos: ${coordStr}</div>
                                <div style="font-size: 0.85em; margin-top: 4px; font-family: monospace;">Op: ${opStr}</div>
                                <div style="margin-top: 4px; font-weight: bold;">${relText}</div>
                              `;
                        grid.appendChild(card);
                    }
                    
                    orderingDiv.appendChild(grid);

                    // Display sublattice translation info for Rank 1 (if computed)
                    if (orderingDiv.hasOwnProperty('sublatticeInfo') || orderingDiv.hasOwnProperty('anyIncomparableLocal')) {
                        const infoDiv = document.createElement('div');
                        infoDiv.style.fontSize = '0.9em';
                        infoDiv.style.marginTop = '6px';
                        if (orderingDiv.anyIncomparableLocal) {
                            infoDiv.innerText = 'Sublattice translation: N/A (non-collinear)';
                        } else if (orderingDiv.sublatticeInfo && Array.isArray(orderingDiv.sublatticeInfo)) {
                            const t = orderingDiv.sublatticeInfo.map(v => decimalToFraction(v)).join(', ');
                            infoDiv.innerText = `Sublattice translation t = (${t}) (fractional)`;
                        } else {
                            // no opposite member found
                            infoDiv.innerText = 'Sublattice translation: N/A';
                        }
                        orderingDiv.appendChild(infoDiv);
                    }

                    // Determine overall ordering status from group colors
                    (() => {
                        // Consider only groups with more than one member
                        const nonSingleGroups = groups.filter(g => g.members && g.members.length > 1);
                        let overallStatus = null; // 'ferro', 'antiferro', 'mixed', null
                        if (nonSingleGroups.length === 0) {
                            overallStatus = 'ncl';
                        } else {
                            let anyGreen = false;
                            let anyRed = false;
                            for (let g of nonSingleGroups) {
                                const rels = g.members.map(mm => mm.relation || mm.relationToRef);
                                let gColor = 'neutral';
                                if (rels.some(r => r === 'Incomparable')) gColor = 'neutral';
                                else if (rels.every(r => r === 'Same')) gColor = 'green';
                                else if (rels.some(r => r === 'Opposite')) gColor = 'red';
                                if (gColor === 'green') anyGreen = true;
                                if (gColor === 'red') anyRed = true;
                            }
                            if (anyGreen) overallStatus = 'ferro';
                            else if (anyRed) overallStatus = 'antiferro';
                            else overallStatus = 'mixed';
                        }

                        const summaryDiv = document.createElement('div');
                        summaryDiv.style.fontWeight = '700';
                        summaryDiv.style.marginTop = '8px';
                        if (overallStatus === 'ferro') {
                            summaryDiv.innerText = 'Ferroically Order';
                            summaryDiv.style.color = '#1b5e20';
                        } else if (overallStatus === 'antiferro') {
                            summaryDiv.innerText = 'Anti-Ferroically Order';
                            summaryDiv.style.color = '#b71c1c';
                        } else if (overallStatus === 'mixed') {
                            summaryDiv.innerText = 'Mixed / Non-collinear ordering';
                            summaryDiv.style.color = '#666';
                        } else if (overallStatus === 'ncl') {
                            summaryDiv.innerText = 'NCL order';
                            summaryDiv.style.color = '#666';
                        } else {
                            summaryDiv.innerText = 'No grouped pairings';
                            summaryDiv.style.color = '#666';
                        }
                        orderingDiv.appendChild(summaryDiv);
                    })();

        // Removed automatic preview of the representative atom; user can click a card to show it
                }
                
                div.appendChild(orderingDiv);
            }
            
            // Use pre-calculated componentGroups for display
            const groups = componentGroups;
            
            // Display Number of Independent Components
            const countDiv = document.createElement('div');
            countDiv.style.marginBottom = "10px";
            countDiv.style.fontStyle = "italic";
            countDiv.innerText = `Number of independent components: ${groups.length}`;
            div.appendChild(countDiv);
            
            // Matrix View with Actual Labels
            const compMap = new Array(nComponents).fill(null); 
            const colorMap = new Array(nComponents).fill(null);
            
            // Generate colors for groups
            const groupColors = [];
            for(let i=0; i<groups.length; i++) {
                // Use HSL to generate distinct colors
                // Avoid very light colors (hard to read on white)
                const hue = (i * 137.508) % 360; 
                groupColors.push(`hsl(${hue}, 80%, 35%)`);
            }

            groups.forEach((g, gIdx) => {
                const refMember = g.members[0];
                const refLabel = labels[refMember.idx];
                const refFactor = refMember.factor;
                
                // Only color if there are relations (more than 1 member in the group)
                const color = (g.members.length > 1) ? groupColors[gIdx] : "black";

                g.members.forEach(m => {
                    // val(m) = (m.factor / refFactor) * val(ref)
                    const ratio = m.factor / refFactor;
                    let s = "";
                    if (Math.abs(ratio - 1) < GROUP_TOL) s = refLabel;
                    else if (Math.abs(ratio + 1) < GROUP_TOL) s = "-" + refLabel;
                    else {
                        const num = ratio.toFixed(3).replace(/\.?0+$/, '');
                        s = num + refLabel;
                    }
                    compMap[m.idx] = s;
                    colorMap[m.idx] = color;
                });
            });
            
            for(let i=0; i<nComponents; i++) {
                if (!compMap[i]) {
                    compMap[i] = "0";
                    colorMap[i] = "#999"; // Grey for zero
                }
            }
            
            // If this is a System (whole-group) result, convert basis to Cartesian for display
            // so the system response is shown in Cartesian coordinates as requested.
            if (!wyckoff && group) {
                try {
                    const latM = getIdealLatticeMatrix(group);
                    // helper: transform tensor components by a 3x3 matrix M (applies M on each index)
                    const transformTensorToCartesianGlobal = (vec, rank, M) => {
                        if (rank === 0) return vec;
                        let out = new Float64Array(vec.length);
                        const n = Math.pow(3, rank);
                        for (let flatOut = 0; flatOut < n; flatOut++) {
                            let idxsOut = [];
                            let rem = flatOut;
                            for (let d = 0; d < rank; d++) {
                                idxsOut.unshift(rem % 3);
                                rem = Math.floor(rem / 3);
                            }
                            let sum = 0;
                            for (let flatIn = 0; flatIn < n; flatIn++) {
                                let idxsIn = [];
                                let remIn = flatIn;
                                for (let d = 0; d < rank; d++) {
                                    idxsIn.unshift(remIn % 3);
                                    remIn = Math.floor(remIn / 3);
                                }
                                let prod = 1;
                                for (let d = 0; d < rank; d++) {
                                    prod *= M[idxsOut[d]][idxsIn[d]];
                                }
                                sum += prod * vec[flatIn];
                            }
                            out[flatOut] = sum;
                        }
                        return out;
                    };

                    // Convert each basis vector (which is in fractional/crystal basis) to Cartesian for display
                    // Build a converted basis array to use when formatting matrices
                    const latM_use = latM;
                    for (let b = 0; b < basis.length; b++) {
                        try {
                            const cart = transformTensorToCartesianGlobal(basis[b], rank, latM_use);
                            // Replace the representative labels mapping values with Cartesian values
                            // We will display using the Cartesian components by swapping compMap source later
                            // To keep minimal changes, write these Cartesian values temporarily into a hidden array
                            basis[b].__cart = cart;
                        } catch(e) {
                            console.warn('Failed to convert basis to Cartesian for system display', e);
                        }
                    }
                } catch(e) {
                    console.warn('Could not compute lattice matrix for Cartesian conversion', e);
                }
            }

            let matHtml = "";
            const useVoigt = document.getElementById('voigtCheck').checked;
            
            if (useVoigt) {
                if (rank === 1) matHtml = formatSymbolicMatrix1(compMap, colorMap); // Rank 1 is same
                else if (rank === 2) matHtml = formatVoigtMatrix2(compMap, colorMap);
                else if (rank === 3) matHtml = formatVoigtMatrix3(compMap, colorMap);
                else if (rank === 4) matHtml = formatVoigtMatrix4(compMap, colorMap);
                else if (rank === 5) matHtml = formatVoigtMatrix5(compMap, colorMap);
            } else {
                if (rank === 1) matHtml = formatSymbolicMatrix1(compMap, colorMap);
                else if (rank === 2) matHtml = formatSymbolicMatrix2(compMap, colorMap);
                else if (rank === 3) matHtml = formatSymbolicMatrix3(compMap, colorMap);
                else if (rank === 4) matHtml = formatSymbolicMatrix4(compMap, colorMap);
                else if (rank === 5) matHtml = formatSymbolicMatrix5(compMap, colorMap);
            }
            
      // The general tensor matrix is shown only for the "System" (whole-group) view.
      if (matHtml && !wyckoff) {
        const matDiv = document.createElement('div');
        matDiv.className = 'matrix-display';
        matDiv.innerHTML = matHtml;
        div.appendChild(matDiv);
      }
            
            container.appendChild(div);
        }

        function formatSymbolicMatrix1(comps, colors) {
            let s = "";
            for(let i=0; i<3; i++) {
                const val = comps[i].padStart(8);
                const col = colors[i];
                s += `| <span style="color:${col}">${val}</span> ; |\n`;
            }
            return s;
        }

        function formatSymbolicMatrix2(comps, colors) {
            let s = "";
            for(let i=0; i<3; i++) {
                s += "| ";
                for(let j=0; j<3; j++) {
                    const idx = i*3+j;
                    const val = comps[idx].padStart(8);
                    const col = colors[idx];
                    s += `<span style="color:${col}">${val}</span>`;
                    if (j < 2) s += "; ";
                }
                s += " |\n";
            }
            return s;
        }
        function exportOrbitJSON() {
            const wSel = document.getElementById('wyckoffSelect').value;

            // reject empty or "whole"
            if (!currentGroup || !currentGroup.wyckoff || wSel === "" || wSel === "whole") {
                alert("Select a Wyckoff position (not 'System') first.");
                return;
            }

            const wyckoff = currentGroup.wyckoff[parseInt(wSel, 10)];
            const orbit = getWyckoffOrbit(currentGroup, wyckoff, wyckoffRandomVars);

            // optional redraw (keep if you want)
            renderOrbitSpheres(orbit);

            // ===== EXPAND ORBIT WITH PERIODIC IMAGES (treat ghosts as real) =====
            const EPS = 1e-3;
            const near = (a, b, eps) => Math.abs(a - b) < eps;

            // EXACT same logic as renderOrbitSpheres()
            const axisShifts = (u) => {
                const s = [0];
                if (near(u, 0, EPS)) s.push(+1);
                if (near(u, 1, EPS)) s.push(-1);
                return s;
            };

            const atomsExpanded = [];
            let id = 1;

            for (let i = 0; i < orbit.length; i++) {
                const o = orbit[i];
                const p = o.coord; // fractional [0..1)

                const sx = axisShifts(p[0]);
                const sy = axisShifts(p[1]);
                const sz = axisShifts(p[2]);

                for (const dx of sx) for (const dy of sy) for (const dz of sz) {
                const pUnwrapped = [p[0] + dx, p[1] + dy, p[2] + dz];

                // keep also wrapped into [0,1)
                const pWrapped = pUnwrapped.map(v => {
                    let x = v - Math.floor(v + 1e-6);
                    if (Math.abs(x - 1.0) < 1e-6) x = 0.0;
                    return x;
                });

                atomsExpanded.push({
                    atom: id++,
                    parent_atom: i + 1,
                    shift: [dx, dy, dz],
                    frac_unwrapped: pUnwrapped,
                    frac: pWrapped,
                    sym: o.coordSym,
                    op: o.op?.str || ""
                });
                }
            }

            // ===== BUILD EXPORT OBJECT =====
            const out = {
                group_index: parseInt(document.getElementById('groupSelect').value, 10),
                group_name: currentGroup.name,
                wyckoff_label: wyckoff.label,
                wyckoff_coord: wyckoff.coord,
                random_vars_used: wyckoffRandomVars,
                atoms: atomsExpanded
            };

            const blob = new Blob([JSON.stringify(out, null, 2)], { type: "application/json" });
            const url = URL.createObjectURL(blob);

            const a = document.createElement("a");
            a.href = url;
            a.download = `orbit_${wyckoff.label.replace(/\s+/g,'_')}.json`;
            document.body.appendChild(a);
            a.click();
            a.remove();

            URL.revokeObjectURL(url);
            }

        function formatSymbolicMatrix3(comps, colors) {
            let s = "";
            const axes = ['x', 'y', 'z'];
            for(let i=0; i<3; i++) {
                s += `i = ${axes[i]}:\n`;
                for(let j=0; j<3; j++) {
                    s += "| ";
                    for(let k=0; k<3; k++) {
                        const idx = i*9 + j*3 + k;
                        const val = comps[idx].padStart(8);
                        const col = colors[idx];
                        s += `<span style="color:${col}">${val}</span>`;
                        if (k < 2) s += "; ";
                    }
                    s += " |\n";
                }
                s += "\n";
            }
            return s;
        }

        function formatSymbolicMatrix4(comps, colors) {
            let s = "";
            const axes = ['x', 'y', 'z'];
            for(let i=0; i<3; i++) {
                for(let j=0; j<3; j++) {
                    s += `i = ${axes[i]}, j = ${axes[j]}:\n`;
                    for(let k=0; k<3; k++) {
                        s += "| ";
                        for(let l=0; l<3; l++) {
                            const idx = i*27 + j*9 + k*3 + l;
                            const val = comps[idx].padStart(8);
                            const col = colors[idx];
                            s += `<span style="color:${col}">${val}</span>`;
                            if (l < 2) s += "; ";
                        }
                        s += " |\n";
                    }
                    s += "\n";
                }
            }
            return s;
        }

        function formatSymbolicMatrix5(comps, colors) {
            let s = "";
            const axes = ['x', 'y', 'z'];
            for(let i=0; i<3; i++) {
                for(let j=0; j<3; j++) {
                    for(let k=0; k<3; k++) {
                        s += `i = ${axes[i]}, j = ${axes[j]}, k = ${axes[k]}:\n`;
                        for(let l=0; l<3; l++) {
                            s += "| ";
                            for(let m=0; m<3; m++) {
                                const idx = i*81 + j*27 + k*9 + l*3 + m;
                                const val = comps[idx].padStart(8);
                                const col = colors[idx];
                                s += `<span style="color:${col}">${val}</span>`;
                                if (m < 2) s += "; ";
                            }
                            s += " |\n";
                        }
                        s += "\n";
                    }
                }
            }
            return s;
        }

        // --- Voigt Formatting ---
        
        // Voigt Map: 1:xx(0), 2:yy(4), 3:zz(8), 4:yz(5), 5:xz(2), 6:xy(1)
        // Note: Standard Voigt uses 2yz, 2xz, 2xy for strain, but usually just components for physical tensors.
        // We use the indices: 0, 4, 8, 5, 2, 1.
        const VOIGT_INDICES = [0, 4, 8, 5, 2, 1];
        const VOIGT_LABELS = ['1', '2', '3', '4', '5', '6'];

        function formatVoigtMatrix2(comps, colors) {
            // Rank 2: Vector of 6
            let s = "| ";
            for(let v=0; v<6; v++) {
                const idx = VOIGT_INDICES[v];
                const val = comps[idx].padStart(8);
                const col = colors[idx];
                s += `<span style="color:${col}">${val}</span> `;
            }
            s += "|\n";
            // Add labels below
            s += "  ";
            for(let v=0; v<6; v++) s += VOIGT_LABELS[v].padStart(8) + " ";
            s += "\n";
            return s;
        }

        function formatVoigtMatrix3(comps, colors) {
            // Rank 3: T_ijk -> T_i_nu. 3x6 Matrix.
            // i is row, nu is col.
            // Index in flat array: i*9 + VOIGT_INDICES[nu]
            let s = "";
            // Header
            s += "     ";
            for(let v=0; v<6; v++) s += VOIGT_LABELS[v].padStart(8) + " ";
            s += "\n";
            
            const rowLabels = ['1', '2', '3'];
            for(let i=0; i<3; i++) {
                s += rowLabels[i] + " | ";
                for(let v=0; v<6; v++) {
                    const idx = i*9 + VOIGT_INDICES[v];
                    const val = comps[idx].padStart(8);
                    const col = colors[idx];
                    s += `<span style="color:${col}">${val}</span> `;
                }
                s += "|\n";
            }
            return s;
        }

        function formatVoigtMatrix4(comps, colors) {
            // Rank 4: T_ijkl -> T_ij_nu. 3x3x6.
            // Display as 3 blocks of 3x6 (i=1..3)
            let s = "";
            const rowLabels = ['1', '2', '3'];
            
            for(let i=0; i<3; i++) {
                s += `i = ${rowLabels[i]}:\n`;
                s += "     ";
                for(let v=0; v<6; v++) s += VOIGT_LABELS[v].padStart(8) + " ";
                s += "\n";
                
                for(let j=0; j<3; j++) {
                    s += rowLabels[j] + " | ";
                    for(let v=0; v<6; v++) {
                        // Index: i*27 + j*9 + VOIGT_INDICES[v]
                        const idx = i*27 + j*9 + VOIGT_INDICES[v];
                        const val = comps[idx].padStart(8);
                        const col = colors[idx];
                        s += `<span style="color:${col}">${val}</span> `;
                    }
                    s += "|\n";
                }
                s += "\n";
            }
            return s;
        }

        function formatVoigtMatrix5(comps, colors) {
            // Rank 5: T_ijklm -> T_ijk_nu. 3x3x3x6.
            let s = "";
            const rowLabels = ['1', '2', '3'];
            
            for(let i=0; i<3; i++) {
                for(let j=0; j<3; j++) {
                    s += `i = ${rowLabels[i]}, j = ${rowLabels[j]}:\n`;
                    s += "     ";
                    for(let v=0; v<6; v++) s += VOIGT_LABELS[v].padStart(8) + " ";
                    s += "\n";
                    
                    for(let k=0; k<3; k++) {
                        s += rowLabels[k] + " | ";
                        for(let v=0; v<6; v++) {
                            // Index: i*81 + j*27 + k*9 + VOIGT_INDICES[v]
                            const idx = i*81 + j*27 + k*9 + VOIGT_INDICES[v];
                            const val = comps[idx].padStart(8);
                            const col = colors[idx];
                            s += `<span style="color:${col}">${val}</span> `;
                        }
                        s += "|\n";
                    }
                    s += "\n";
                }
            }
            return s;
        }

        function getTensorLabels_local(rank) {
            const axes = ['x', 'y', 'z'];
            if (rank === 1) return axes;
            
            let labels = axes;
            for (let r = 2; r <= rank; r++) {
                const newLabels = [];
                for (let l of labels) {
                    for (let a of axes) {
                        newLabels.push(l + a);
                    }
                }
                labels = newLabels;
            }
            return labels;
        }

        // --- Linear Algebra Utils ---

        function matmul(A, B) {
            const C = [[0,0,0],[0,0,0],[0,0,0]];
            for(let i=0; i<3; i++)
                for(let j=0; j<3; j++)
                    for(let k=0; k<3; k++)
                        C[i][j] += A[i][k] * B[k][j];
            return C;
        }

        function matvec(A, v) {
            const r = [0,0,0];
            for(let i=0; i<3; i++)
                for(let j=0; j<3; j++)
                    r[i] += A[i][j] * v[j];
            return r;
        }

        function vecAdd(a, b) { return [a[0]+b[0], a[1]+b[1], a[2]+b[2]]; }
        function vecSub(a, b) { return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]; }
        
        function vecDot(a, b) {
            let s = 0;
            for(let i=0; i<a.length; i++) s += a[i]*b[i];
            return s;
        }
        
        function vecComb(c1, v1, c2, v2) {
            const r = new Float64Array(v1.length);
            for(let i=0; i<v1.length; i++) r[i] = c1*v1[i] + c2*v2[i];
            return r;
        }

        function matrixEqual(A, B) {
            for(let i=0; i<3; i++)
                for(let j=0; j<3; j++)
                    if (A[i][j] !== B[i][j]) return false;
            return true;
        }

        function isIntegerVec(v) {
            return Math.abs(v[0] - Math.round(v[0])) < 1e-5 &&
                   Math.abs(v[1] - Math.round(v[1])) < 1e-5 &&
                   Math.abs(v[2] - Math.round(v[2])) < 1e-5;
        }

        function det3x3(m) {
            return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                   m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                   m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        }

        function getIdealLatticeMatrix(group) {
              if (CUSTOM_CELL_GLOBAL_ENABLED && CUSTOM_CELL_GLOBAL) {
                const basis = latticeBasisFromParams(CUSTOM_CELL_GLOBAL);
                if (basis) {
                  return [
                    [basis.a[0], basis.b[0], basis.c[0]],
                    [basis.a[1], basis.b[1], basis.c[1]],
                    [basis.a[2], basis.b[2], basis.c[2]]
                  ];
                }
              }
            // Determine system from group number
            const name = group.name; 
            const match = name.match(/BNS:\s*(\d+)\./);
            let num = 1;
            if (match) num = parseInt(match[1]);

            const system = crystalSystemFromSGNumber(num);
            if ((system === "triclinic" || system === "monoclinic") &&
                CUSTOM_CELL_ENABLED && CUSTOM_CELL_ENABLED[system]) {
                const custom = CUSTOM_CELL_PARAMS && CUSTOM_CELL_PARAMS[system];
                const basis = custom ? latticeBasisFromParams(custom) : null;
                if (basis) {
                    return [
                        [basis.a[0], basis.b[0], basis.c[0]],
                        [basis.a[1], basis.b[1], basis.c[1]],
                        [basis.a[2], basis.b[2], basis.c[2]]
                    ];
                }
            }
            
            let a=1, b=1, c=1;
            let alpha=90, beta=90, gamma=90;

            if (num >= 195) { // Cubic
                // a=b=c, 90,90,90
            } else if (num >= 168) { // Hexagonal
                // a=b, 90,90,120
                gamma = 120;
            } else if (num >= 143) { // Trigonal
                // Usually Hex setting in BNS/OG for these ranges?
                // Assume Hex setting for standard
                gamma = 120;
            } else if (num >= 75) { // Tetragonal
                // a=b, 90,90,90
            } else if (num >= 16) { // Orthorhombic
                // a!=b!=c, 90,90,90
                // Ratios don't matter for point group
            } else if (num >= 3) { // Monoclinic
                // 90, beta, 90. 
                // Beta doesn't affect point group matrices in standard setting
            } else { // Triclinic
                // All free.
            }

            // Convert to radians
            const ar = alpha * Math.PI / 180;
            const br = beta * Math.PI / 180;
            const gr = gamma * Math.PI / 180;

            const ca = Math.cos(ar);
            const cb = Math.cos(br);
            const cg = Math.cos(gr);
            const sg = Math.sin(gr);
            
            const V_term = Math.sqrt(Math.max(0, 1 - ca*ca - cb*cb - cg*cg + 2*ca*cb*cg));
            
            const M = [
                [a, b*cg, c*cb],
                [0, b*sg, c*(ca - cb*cg)/sg],
                [0, 0,    c*V_term/sg]
            ];
            return M;
        }

        

        function invertMatrix3x3(m) {
            const det = det3x3(m);
            if (Math.abs(det) < 1e-9) return null;
            const invDet = 1.0 / det;
            
            const inv = [[0,0,0],[0,0,0],[0,0,0]];
            
            inv[0][0] = (m[1][1]*m[2][2] - m[1][2]*m[2][1]) * invDet;
            inv[0][1] = (m[0][2]*m[2][1] - m[0][1]*m[2][2]) * invDet;
            inv[0][2] = (m[0][1]*m[1][2] - m[0][2]*m[1][1]) * invDet;
            
            inv[1][0] = (m[1][2]*m[2][0] - m[1][0]*m[2][2]) * invDet;
            inv[1][1] = (m[0][0]*m[2][2] - m[0][2]*m[2][0]) * invDet;
            inv[1][2] = (m[0][2]*m[1][0] - m[0][0]*m[1][2]) * invDet;
            
            inv[2][0] = (m[1][0]*m[2][1] - m[1][1]*m[2][0]) * invDet;
            inv[2][1] = (m[0][1]*m[2][0] - m[0][0]*m[2][1]) * invDet;
            inv[2][2] = (m[0][0]*m[1][1] - m[0][1]*m[1][0]) * invDet;
            
            return inv;
        }

        // --- 3D Viewer (Three.js) ---
// NOTE: scene3d/camera3d/renderer3d/controls3d/orbitGroup3d are declared near the top of the file.

function initViewer3D() {
    const container = document.getElementById("viewer3d");
    if (!container) return;

    // If already initialized, just ensure size is correct
    if (renderer3d && camera3d) {
        const ww = container.clientWidth;
        const hh = container.clientHeight;
        if (ww && hh) {
            camera3d.aspect = ww / hh;
            camera3d.updateProjectionMatrix();
            renderer3d.setSize(ww, hh);
        }
        return;
    }

    // If the container is not sized yet, retry next frame.
    const w = container.clientWidth;
    const h = container.clientHeight;
    if (!w || !h) {
        requestAnimationFrame(initViewer3D);
        return;
    }

    scene3d = new THREE.Scene();
    scene3d.background = new THREE.Color(0xffffff);
    


    camera3d = new THREE.PerspectiveCamera(45, w / h, 0.01, 200);
    camera3d.position.set(1.6, 1.2, 1.8);
    camera3d.lookAt(0, 0, 0);

    renderer3d = new THREE.WebGLRenderer({ antialias: true, alpha: false });
    renderer3d.setPixelRatio(window.devicePixelRatio || 1);
    renderer3d.setSize(w, h);

    // Mount canvas (clear first so we don't stack canvases)
    container.innerHTML = "";
    container.appendChild(renderer3d.domElement);
    // --- tooltip init + hover events (bind once) ---
    if (!atomTooltipEl) atomTooltipEl = document.getElementById("atom-tooltip");
    if (!atomRaycaster) atomRaycaster = new THREE.Raycaster();
    if (!atomMouseNdc) atomMouseNdc = new THREE.Vector2();

    if (!__atomHoverBound) {
    __atomHoverBound = true;

    renderer3d.domElement.addEventListener("mousemove", onAtomHoverMove, false);
    renderer3d.domElement.addEventListener("mouseleave", () => {
        if (atomTooltipEl) atomTooltipEl.style.display = "none";
    }, false);
    }


    // Prevent browser gestures from hijacking mouse/touch (ant for trackpads)
    renderer3d.domElement.style.touchAction = "none";
    container.oncontextmenu = (e) => e.preventDefault();

    // lights
    scene3d.add(new THREE.AmbientLight(0xffffff, 0.55));
    const dir = new THREE.DirectionalLight(0xffffff, 0.9);
    dir.position.set(2, 3, 4);
    scene3d.add(dir);

    // controls
    // controls (support both legacy THREE.OrbitControls and global OrbitControls)
    const OrbitControlsCtor = THREE.OrbitControls || window.OrbitControls;
    addCrystalAxes({ scale: 0.45 });
    // Apply unit cell that was requested before the viewer was ready
    if (__PENDING_UNITCELL_NAME) {
    setUnitCellFromGroupName(__PENDING_UNITCELL_NAME);
    }


    if (OrbitControlsCtor) {
        controls3d = new OrbitControlsCtor(camera3d, renderer3d.domElement);
        controls3d.enableDamping = true;
        controls3d.dampingFactor = 0.08;

        controls3d.enableRotate = true;
        controls3d.enableZoom = true;
        controls3d.enablePan = true;

        controls3d.mouseButtons = {
            LEFT: THREE.MOUSE.ROTATE,
            MIDDLE: THREE.MOUSE.DOLLY,
            RIGHT: THREE.MOUSE.PAN
        };

        controls3d.zoomSpeed = 1.0;
        controls3d.panSpeed = 1.0;
        controls3d.target.set(0, 0, 0);
        controls3d.update();
    } else {
        controls3d = null;
        console.warn("OrbitControls not found. OrbitControls.js didn't load or isn't compatible.");
    }

    // ---- dynamic unit cell wireframe (changes with selected space-group family) ----
    

    // Parse "BNS: 63.462 ..." -> 63
    

    // One resize handler only
    if (!window.__viewer3dResizeBound) {
        window.__viewer3dResizeBound = true;
        window.addEventListener("resize", () => {
            const c = document.getElementById("viewer3d");
            if (!c || !renderer3d || !camera3d) return;
            const ww = c.clientWidth;
            const hh = c.clientHeight;
            if (!ww || !hh) return;
            camera3d.aspect = ww / hh;
            camera3d.updateProjectionMatrix();
            renderer3d.setSize(ww, hh);
        });
    }

    // animation loop (start once)
    if (!window.__viewer3dAnimRunning) {
        window.__viewer3dAnimRunning = true;
        (function animate() {
            requestAnimationFrame(animate);
            if (controls3d) controls3d.update();
            if (renderer3d && scene3d && camera3d) renderer3d.render(scene3d, camera3d);
        })();
    }
}
// ==================== UNIT CELL (MAILLE) LOGIC ====================
let cellWireframe3d = null;
let __PENDING_UNITCELL_NAME = null;

// Parse "BNS: 63.462 ..." -> 63
function parseBNSNumberMajor(groupName) {
  const m = String(groupName).match(/BNS:\s*([0-9]+)\./);
  return m ? parseInt(m[1], 10) : null;
}
function addXYZAxes(size = 0.8) {
  if (!scene3d) return;

  // remove old axes if any
  if (axesHelper3d) {
    scene3d.remove(axesHelper3d);
    axesHelper3d = null;
  }

  axesHelper3d = new THREE.AxesHelper(size);
  scene3d.add(axesHelper3d);
}

// crystal system from International space-group number ranges
function crystalSystemFromSGNumber(n) {
  if (n === null || isNaN(n)) return "cubic";
  if (n >= 1 && n <= 2) return "triclinic";
  if (n >= 3 && n <= 15) return "monoclinic";
  if (n >= 16 && n <= 74) return "orthorhombic";
  if (n >= 75 && n <= 142) return "tetragonal";
  if (n >= 143 && n <= 167) return "trigonal";
  if (n >= 168 && n <= 194) return "hexagonal";
  return "cubic";
}
function latticeBasisFromParams(params) {
  const deg = Math.PI / 180;
  const a = params.a;
  const b = params.b;
  const c = params.c;
  const alpha = params.alpha * deg;
  const beta = params.beta * deg;
  const gamma = params.gamma * deg;

  const cosAlpha = Math.cos(alpha);
  const cosBeta = Math.cos(beta);
  const cosGamma = Math.cos(gamma);
  const sinGamma = Math.sin(gamma);

  if (Math.abs(sinGamma) < 1e-6) return null;

  const aVec = [a, 0, 0];
  const bVec = [b * cosGamma, b * sinGamma, 0];
  const cX = c * cosBeta;
  const cY = c * (cosAlpha - cosBeta * cosGamma) / sinGamma;
  const cZSq = Math.max(c * c - cX * cX - cY * cY, 0);
  const cZ = Math.sqrt(cZSq);

  const cVec = [cX, cY, cZ];

  const cross = [
    bVec[1] * cVec[2] - bVec[2] * cVec[1],
    bVec[2] * cVec[0] - bVec[0] * cVec[2],
    bVec[0] * cVec[1] - bVec[1] * cVec[0]
  ];
  const volume = Math.abs(aVec[0] * cross[0] + aVec[1] * cross[1] + aVec[2] * cross[2]);
  if (volume > 0) {
    const scale = Math.pow(1 / volume, 1 / 3);
    return {
      a: aVec.map(v => v * scale),
      b: bVec.map(v => v * scale),
      c: cVec.map(v => v * scale)
    };
  }

  return { a: aVec, b: bVec, c: cVec };
}
function getCellBasis(system) {
  let a = [1, 0, 0], b = [0, 1, 0], c = [0, 0, 1];

  if (CUSTOM_CELL_GLOBAL_ENABLED && CUSTOM_CELL_GLOBAL) {
    const basis = latticeBasisFromParams(CUSTOM_CELL_GLOBAL);
    if (basis) return basis;
  }

  if ((system === "triclinic" || system === "monoclinic") &&
      CUSTOM_CELL_ENABLED && CUSTOM_CELL_ENABLED[system]) {
    const custom = CUSTOM_CELL_PARAMS && CUSTOM_CELL_PARAMS[system];
    if (custom) {
      const basis = latticeBasisFromParams(custom);
      if (basis) return basis;
    }
  }

  if (system === "tetragonal") {
    c = [0, 0, 1.4];
  } else if (system === "orthorhombic") {
    a = [1.2, 0, 0];
    b = [0, 0.8, 0];
    c = [0, 0, 1.4];
  } else if (system === "hexagonal" || system === "trigonal") {
    a = [1, 0, 0];
    b = [-0.5, Math.sqrt(3)/2, 0];
    c = [0, 0, 1.4];
  } else if (system === "monoclinic") {
    a = [1.2, 0, 0];
    b = [0, 1.0, 0];
    c = [0.35, 0, 1.3];
  } else if (system === "triclinic") {
    a = [1.1, 0.1, 0];
    b = [0.2, 0.95, 0.1];
    c = [0.3, 0.2, 1.1];
  }
  return { a, b, c };
}

function buildCellWireframe(system) {


  const V = [];
  const addEdge = (a, b) => { V.push(...a, ...b); };

  // Use the shared helper to ensure consistency
  let { a, b, c } = getCellBasis(system);

  // (Removed duplicated if/else logic here)

  const O  = [0, 0, 0];
  const A  = a;
  const B  = b;
  const AB = [a[0]+b[0], a[1]+b[1], a[2]+b[2]];
  const C  = c;
  const AC = [a[0]+c[0], a[1]+c[1], a[2]+c[2]];
  const BC = [b[0]+c[0], b[1]+c[1], b[2]+c[2]];
  const ABC= [a[0]+b[0]+c[0], a[1]+b[1]+c[1], a[2]+b[2]+c[2]];

  const center = [ABC[0]/2, ABC[1]/2, ABC[2]/2];
  const shift = (p) => [p[0]-center[0], p[1]-center[1], p[2]-center[2]];

  const o  = shift(O),  aa = shift(A),  bb = shift(B),  ab = shift(AB);
  const cc = shift(C),  ac = shift(AC), bc = shift(BC), abc = shift(ABC);

  addEdge(o, aa); addEdge(o, bb); addEdge(o, cc);
  addEdge(aa, ab); addEdge(aa, ac);
  addEdge(bb, ab); addEdge(bb, bc);
  addEdge(cc, ac); addEdge(cc, bc);
  addEdge(ab, abc); addEdge(ac, abc); addEdge(bc, abc);

  const geo = new THREE.BufferGeometry();
  geo.setAttribute("position", new THREE.Float32BufferAttribute(V, 3));
  const mat = new THREE.LineBasicMaterial({ color: 0x000000 });
  CELL_BASIS = { a, b, c };
  return new THREE.LineSegments(geo, mat);
}
// Build P such that: W(harmonics) = P * c(basis coeffs)
// P has shape (nHarm x nBasis)
function buildBasisToHarmonicMatrix(fullRank, basis, spinIndex, harmonics) {
  const nH = harmonics.length;
  const nB = basis.length;
  const P = Array.from({length:nH}, () => Array(nB).fill(0));

  for (let b = 0; b < nB; b++) {
    const { spatialRank, spatialTensor } = getSpatialTensorForProjection(
      fullRank,
      spinIndex,
      basis[b],
      (CELL_BASIS && CELL_BASIS.a) ? [
        [CELL_BASIS.a[0], CELL_BASIS.b[0], CELL_BASIS.c[0]],
        [CELL_BASIS.a[1], CELL_BASIS.b[1], CELL_BASIS.c[1]],
        [CELL_BASIS.a[2], CELL_BASIS.b[2], CELL_BASIS.c[2]]
      ] : [[1,0,0],[0,1,0],[0,0,1]]
    );
    const w = projectTensorToSpecificHarmonics(
      spatialRank,
      new Float64Array(spatialTensor),
      harmonics,
      { Nth: PROJ_NTH, Nph: PROJ_NPH }
    );
    for (let h = 0; h < nH; h++) P[h][b] = w[h];
  }
  return P;
}

function clearCellWireframe() {
  if (cellWireframe3d && scene3d) {
    scene3d.remove(cellWireframe3d);
    cellWireframe3d.geometry.dispose();
    cellWireframe3d.material.dispose();
    cellWireframe3d = null;
  }
}

function setUnitCellFromGroupName(groupName) {
  // viewer not ready yet -> remember request
  const major = parseBNSNumberMajor(groupName);
  const system = crystalSystemFromSGNumber(major);
  if (!scene3d) {
    __PENDING_UNITCELL_NAME = groupName;
    maybeOpenCellParamsModal(system);
    return;
  }
  __PENDING_UNITCELL_NAME = null;

  clearCellWireframe();
  cellWireframe3d = buildCellWireframe(system);
  scene3d.add(cellWireframe3d);
  addCrystalAxes({ scale: 0.45 });
  maybeOpenCellParamsModal(system);
  if (_lastOrbitForTensor) {
    renderTensorIsosurfacesFromSiteTensor(_lastOrbitForTensor);
  }
}
function resetTensorStateForNewPosition() {
  resetTensorSliderMemory();   // clears TV_STATE + ALLOWED_HARM_CACHE
  clearTensorSlidersUI();      // clears the slider DOM

  // IMPORTANT for Bug 8 logic: also clear persisted variables
  for (const k in GLOBAL_BASIS_COEFFS) delete GLOBAL_BASIS_COEFFS[k];
  // Clear the global variable registry too so new Wyckoff selections start fresh
  for (const k in VARIABLE_REGISTRY) delete VARIABLE_REGISTRY[k];
}

function resetAllForWyckoffChange() {
  // Clear computed results
  lastResults = [];
  window.lastResults = lastResults;

  const resultsDiv = document.getElementById('results');
  if (resultsDiv) resultsDiv.innerHTML = '';

  const status = document.getElementById('status');
  if (status) status.innerText = '';

  // Reset tensor state/caches
  resetTensorStateForNewPosition();

  // Clear 3D viewer state
  clearViewerOrbit();
  clearIsoGroup();
  _lastOrbitForTensor = null;
}

function resetTensorSliderMemory() {
  TV_STATE = {};
  for (const k in ALLOWED_HARM_CACHE) delete ALLOWED_HARM_CACHE[k];
}

function clearTensorSlidersUI() {
  const box = document.getElementById('tensor-sliders');
  if (box) box.innerHTML = '';

  const forbidden = document.getElementById('tensor-forbidden');
  if (forbidden) {
    forbidden.style.display = 'none';
    forbidden.innerText = '';
  }

  // ❌ remove this line:
  // TV_STATE = {};
}


function clearIsoGroup() {
  if (!scene3d) return;
  if (isoGroup3d) {
    scene3d.remove(isoGroup3d);
    isoGroup3d = null;
  }
}
function onAtomHoverMove(e) {
  if (!renderer3d || !camera3d || !orbitGroup3d || !atomTooltipEl) return;

  const rect = renderer3d.domElement.getBoundingClientRect();
  const x = ((e.clientX - rect.left) / rect.width) * 2 - 1;
  const y = -(((e.clientY - rect.top) / rect.height) * 2 - 1);

  atomMouseNdc.set(x, y);
  atomRaycaster.setFromCamera(atomMouseNdc, camera3d);

  // only check spheres currently in orbitGroup3d
  const hits = atomRaycaster.intersectObjects(orbitGroup3d.children, true);

  const hit = hits.find(h => h.object && h.object.userData && h.object.userData.atomTooltip);
  if (!hit) {
    atomTooltipEl.style.display = "none";
    return;
  }

  atomTooltipEl.textContent = hit.object.userData.atomTooltip;

  atomTooltipEl.style.left = (e.clientX + 12) + "px";
  atomTooltipEl.style.top  = (e.clientY + 12) + "px";
  atomTooltipEl.style.display = "block";
}

// Update tooltips AFTER ordering is known (UP/DOWN/NCL)
function updateOrbitSphereTooltips(orbit, groups) {
  if (!orbitGroup3d || !Array.isArray(orbit)) return;

  // relation per atomIdx (0-based)
  const rel = new Array(orbit.length).fill("NCL");
  if (Array.isArray(groups)) {
    groups.forEach(g => {
      (g.members || []).forEach(m => {
        const idx = m.atomIdx;
        const r = m.relationToRef || m.relation;
        if (idx != null && r) rel[idx] = r;
      });
    });
  }

  orbitGroup3d.traverse(obj => {
    if (!obj.isMesh) return;
    const i = obj.userData?.atomBaseIdx;
    if (i == null) return;

    const p = orbit[i];
    if (!p) return;

    const atomLabel = `Atom ${i + 1}`;
    const pos = p.coordSym || p.coord || "(unknown)";
    const op  = p.op?.str || p.opStr || "(unknown)";
    const relText = (rel[i] === "Same") ? "UP (+)" : (rel[i] === "Opposite") ? "DOWN (-)" : "NCL";

    obj.userData.atomTooltip =
      `${atomLabel}\n` +
      `Pos: ${pos}\n` +
      `Op: ${op}\n` +
      `${relText}`;
  });
}

function getSelectedIsoRanks() {
  const el = document.getElementById("rankIsoSelect");
  if (!el) return [];
  const val = parseInt(el.value, 10);
  return Number.isNaN(val) ? [] : [val];
}

function clearViewerOrbit() {
    if (!scene3d) return;
    if (orbitGroup3d) {
        scene3d.remove(orbitGroup3d);
        orbitGroup3d = null;
    }
}
function expandOrbitForViewer(orbit) {
  const EPS = 1e-3;
  const near = (a,b,eps) => Math.abs(a-b) < eps;

  const axisShifts = (u) => {
    const s = [0];
    if (near(u, 0, EPS)) s.push(+1);
    if (near(u, 1, EPS)) s.push(-1);
    return s;
  };

  const out = [];
  for (let i = 0; i < orbit.length; i++) {
    const atom = orbit[i];
    const p = atom.coord; // [0..1)

    const sx = axisShifts(p[0]);
    const sy = axisShifts(p[1]);
    const sz = axisShifts(p[2]);

    for (const dx of sx) for (const dy of sy) for (const dz of sz) {
      const isGhost = (dx !== 0 || dy !== 0 || dz !== 0);

      // IMPORTANT: use the *same* coordinates used by the spheres (p + shift)
      const pu = [p[0] + dx, p[1] + dy, p[2] + dz];

      // same “keep only one unit cell cube” rule as renderOrbitSpheres()
      const x = pu[0] - 0.5;
      const y = pu[1] - 0.5;
      const z = pu[2] - 0.5;
      if (x < -0.5001 || x > 0.5001) continue;
      if (y < -0.5001 || y > 0.5001) continue;
      if (z < -0.5001 || z > 0.5001) continue;

      out.push({
        ...atom,
        originalIndex: i,
        coord: pu,                 // unwrapped (matches sphere placement)
        shift: [dx, dy, dz],
        isGhost
      });
    }
  }
  return out;
}
function matT(A) {
  // A: rows x cols -> cols x rows
  const r = A.length, c = A[0].length;
  const out = Array.from({length:c}, ()=>Array(r).fill(0));
  for (let i=0;i<r;i++) for (let j=0;j<c;j++) out[j][i]=A[i][j];
  return out;
}

function matMul(A,B) {
  const r=A.length, k=A[0].length, c=B[0].length;
  const out = Array.from({length:r}, ()=>Array(c).fill(0));
  for (let i=0;i<r;i++) for (let j=0;j<c;j++) {
    let s=0;
    for (let t=0;t<k;t++) s += A[i][t]*B[t][j];
    out[i][j]=s;
  }
  return out;
}

function matVec(A, v) {
  const r=A.length, c=A[0].length;
  const out = new Array(r).fill(0);
  for (let i=0;i<r;i++) {
    let s=0;
    for (let j=0;j<c;j++) s += A[i][j]*v[j];
    out[i]=s;
  }
  return out;
}

function invertMatrixN(A) {
  const n=A.length;
  // augment with identity
  const M = A.map((row,i)=> row.slice().concat(Array.from({length:n},(_,j)=>i===j?1:0)));
  for (let i=0;i<n;i++) {
    // pivot
    let piv = M[i][i];
    if (Math.abs(piv) < 1e-12) {
      let swap=i+1;
      while (swap<n && Math.abs(M[swap][i])<1e-12) swap++;
      if (swap===n) throw new Error("Singular matrix");
      [M[i],M[swap]]=[M[swap],M[i]];
      piv = M[i][i];
    }
    // normalize row
    for (let j=0;j<2*n;j++) M[i][j] /= piv;
    // eliminate others
    for (let k=0;k<n;k++) if (k!==i) {
      const f = M[k][i];
      for (let j=0;j<2*n;j++) M[k][j] -= f*M[i][j];
    }
  }
  return M.map(row=>row.slice(n));
}

// Least squares solve x = argmin ||A x - b||
// Using normal equations with tiny ridge for stability.
function leastSquares(A, b, ridge = 1e-8) {
  const At = matT(A);
  const AtA = matMul(At, A);
  for (let i=0;i<AtA.length;i++) AtA[i][i] += ridge;
  const inv = invertMatrixN(AtA);
  const Atb = matVec(At, b);
  return matVec(inv, Atb);
}

function fracToCart(u, v, w) {
  const { a, b, c } = CELL_BASIS || { a:[1,0,0], b:[0,1,0], c:[0,0,1] };
  const uu = u - 0.5, vv = v - 0.5, ww = w - 0.5;
  return [
    uu*a[0] + vv*b[0] + ww*c[0],
    uu*a[1] + vv*b[1] + ww*c[1],
    uu*a[2] + vv*b[2] + ww*c[2],
  ];
}

// Convert fractional rotation R (acting on frac coords) to cart rotation:
// R_cart = B * R_frac * B^{-1}, where B=[a b c] (columns)
function fracRotToCartMatrix4(R) {
  const { a, b, c } = CELL_BASIS || { a:[1,0,0], b:[0,1,0], c:[0,0,1] };

  const B = new THREE.Matrix3().set(
    a[0], b[0], c[0],
    a[1], b[1], c[1],
    a[2], b[2], c[2],
  );

  const Binv = B.clone();
  if (Math.abs(Binv.determinant()) < 1e-12) return new THREE.Matrix4(); // fallback identity
  Binv.invert();

  const Rf = new THREE.Matrix3().set(
    R[0][0], R[0][1], R[0][2],
    R[1][0], R[1][1], R[1][2],
    R[2][0], R[2][1], R[2][2],
  );

  const Rc = B.clone().multiply(Rf).multiply(Binv);

  const m4 = new THREE.Matrix4();
  m4.set(
    Rc.elements[0], Rc.elements[1], Rc.elements[2], 0,
    Rc.elements[3], Rc.elements[4], Rc.elements[5], 0,
    Rc.elements[6], Rc.elements[7], Rc.elements[8], 0,
    0, 0, 0, 1
  );
  return m4;
}

// Transform a spatial tensor (rank = spatialRank, flattened) from fractional
// coordinates to Cartesian coordinates using 3x3 matrix M where M[i][j]
// maps fractional index j -> Cartesian index i: out_i = sum_j M[i][j] * in_j
function transformSpatialTensorToCartesian(tensorFlat, rank, M) {
  if (rank === 0) return tensorFlat.slice();
  const n = Math.pow(3, rank);
  const out = new Float64Array(n);

  // iterate over all output indices
  for (let flatOut = 0; flatOut < n; flatOut++) {
    // decode index digits
    let idxsOut = [];
    let tmp = flatOut;
    for (let d = 0; d < rank; d++) { idxsOut.unshift(tmp % 3); tmp = Math.floor(tmp / 3); }

    let sum = 0;
    // sum over all input combinations
    for (let flatIn = 0; flatIn < n; flatIn++) {
      let idxsIn = [];
      let t = flatIn;
      for (let d = 0; d < rank; d++) { idxsIn.unshift(t % 3); t = Math.floor(t / 3); }

      let prod = 1;
      for (let d = 0; d < rank; d++) {
        prod *= M[ idxsOut[d] ][ idxsIn[d] ];
      }
      sum += prod * tensorFlat[flatIn];
    }
    out[flatOut] = sum;
  }

  return out;
}

// Derive the spatial tensor used for harmonic projection.
// For spatialRank=0, compute the cartesian component (Mx/My/Mz) from the full rank-1 tensor.
function getSpatialTensorForProjection(fullRank, spin, tensorFlat, M) {
  const mode = document.getElementById('multipoleTypeSelect')?.value || 'magnetic';

  // Generic handling: Transform the full tensor from Fractional to Cartesian first.
  // This ensures that when we slice by 'spin', we are selecting Cartesian components (x,y,z)
  // rather than Fractional components (a,b,c).
  // This is crucial for non-cubic systems where fractional components do not correspond to orthogonal axes.
  const cartTensor = transformSpatialTensorToCartesian(tensorFlat, fullRank, M);

  if (mode !== 'magnetic') {
    return { spatialRank: fullRank, spatialTensor: cartTensor };
  }

  // After transforming, the first index (spin) corresponds to Cartesian X, Y, Z.
  const spatialCart = sliceSpinComponent(cartTensor, fullRank, spin);
  
  const spatialRank = Math.max(0, fullRank - 1);
  return { spatialRank, spatialTensor: spatialCart };
}

function near(a, b, eps) { return Math.abs(a - b) < eps; }
function renderOrbitSpheres(orbit) {
  initViewer3D();
  if (!scene3d) return;

  clearViewerOrbit();
  orbitGroup3d = new THREE.Group();

  const r = 0.04;
  const sphereGeo = new THREE.SphereGeometry(r, 24, 24);
  const EPS = 1e-3;

  // Convert fractional -> cartesian using current cell basis
  // Requires: CELL_BASIS = { a:[...], b:[...], c:[...] } set elsewhere (same as wireframe)
  function fracToCart(u, v, w) {
    const { a, b, c } = CELL_BASIS || { a:[1,0,0], b:[0,1,0], c:[0,0,1] };

    // center in the cell (same visual convention as your wireframe centered at origin)
    const uu = u - 0.5;
    const vv = v - 0.5;
    const ww = w - 0.5;

    return [
      uu*a[0] + vv*b[0] + ww*c[0],
      uu*a[1] + vv*b[1] + ww*c[1],
      uu*a[2] + vv*b[2] + ww*c[2],
    ];
  }
  function wrap01(u) {
    u = u % 1;
    return (u < 0) ? (u + 1) : u;
    }

  // Decide which periodic images to draw along one axis
  function axisShifts(u) {
    const s = [0];
    if (near(u, 0, EPS)) s.push(+1);
    if (near(u, 1, EPS)) s.push(-1);
    return s;
  }

  // Keep ONLY one unit cell (in fractional space)
  function inUnitCellFrac(u) {
    return (u > -EPS && u < 1 + EPS);
  }

  for (let i = 0; i < orbit.length; i++) {
    const p0 = orbit[i].coord; 
    const p = [wrap01(p0[0]), wrap01(p0[1]), wrap01(p0[2])];

     // fractional [0..1)

    const sx = axisShifts(p[0]);
    const sy = axisShifts(p[1]);
    const sz = axisShifts(p[2]);

    const col = new THREE.Color();
    col.setHSL(0.58, 0.9, 0.55);

    const matMain = new THREE.MeshStandardMaterial({
    color: col,
    transparent: true,   // ✅ enable transparency
    opacity: 0.25,       // 👈 adjust (0.15–0.3 looks good)
    roughness: 0.2,
    metalness: 0.0,

    depthWrite: false,   // 🔥 THIS IS THE KEY LINE
    side: THREE.DoubleSide
    });

    // ghost atoms: same material (or you can lower opacity more if you want)
    const matGhost = matMain;



    for (const dx of sx) for (const dy of sy) for (const dz of sz) {
      const isGhost = (dx !== 0 || dy !== 0 || dz !== 0);

      // shifted fractional coords
      const pu = p[0] + dx;
      const pv = p[1] + dy;
      const pw = p[2] + dz;

      // Keep ONLY one unit cell (fractional check, not cartesian cube check)
      if (!inUnitCellFrac(pu) || !inUnitCellFrac(pv) || !inUnitCellFrac(pw)) continue;

      // Convert to cartesian consistent with the wireframe cell
      const [x, y, z] = fracToCart(pu, pv, pw);

      const sphere = new THREE.Mesh(
      sphereGeo,
      isGhost ? matGhost : matMain
      );
      sphere.userData.atomBaseIdx = i; // <-- IMPORTANT: i is the orbit index in your loop

    // initial tooltip (relation may be updated later)
     const atomLabel = `Atom ${i + 1}`;
     const pos = orbit[i]?.coordSym || orbit[i]?.coord || "(unknown)";
     const op  = orbit[i]?.op?.str || "(unknown)";
     sphere.userData.atomTooltip =
     `${atomLabel}\n` +
     `Pos: ${pos}\n` +
     `Op: ${op}\n` +
     `NCL`;


      sphere.renderOrder = -1;   // ✅ always behind tensors
      sphere.position.set(x, y, z);
      orbitGroup3d.add(sphere);

    }
  }

  scene3d.add(orbitGroup3d);

  if (controls3d) {
    controls3d.target.set(0, 0, 0);
    controls3d.update();
  }
}

// ==================== Ylm math (FROM tensor_visualizer.html) ====================

function factorial(n) {
  if (n <= 1) return 1;
  return n * factorial(n - 1);
}

function legendreP(l, m, x) {
  if (m < 0) {
    const sign = (Math.abs(m) % 2 === 1) ? -1 : 1;
    return sign * factorial(l + m) / factorial(l - Math.abs(m)) * legendreP(l, -m, x);
  }

  let pmm = 1.0;
  if (m > 0) {
    const somx2 = Math.sqrt((1.0 - x) * (1.0 + x));
    let fact = 1.0;
    for (let i = 1; i <= m; i++) {
      pmm *= -fact * somx2;
      fact += 2.0;
    }
  }
  if (l === m) return pmm;

  let pmmp1 = x * (2.0 * m + 1.0) * pmm;
  if (l === m + 1) return pmmp1;

  let pll = 0.0;
  for (let ll = m + 2; ll <= l; ll++) {
    pll = ((2.0 * ll - 1.0) * x * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m);
    pmm = pmmp1;
    pmmp1 = pll;
  }
  return pll;
}

function getRealYlm(l, m, theta, phi) {
  const cosTheta = Math.cos(theta);
  const P = legendreP(l, Math.abs(m), cosTheta);

  const norm =
    Math.sqrt(
      ((2 * l + 1) * factorial(l - Math.abs(m))) /
      (4 * Math.PI * factorial(l + Math.abs(m)))
    );

  const Y = norm * P;

  if (m === 0) return Y;
  if (m > 0) return Math.sqrt(2) * Y * Math.cos(m * phi);
  return Math.sqrt(2) * Y * Math.sin(Math.abs(m) * phi);
}

// ======================================================================



// --------- local tensor -> Ylm projection (numerical) ---------
function _idxToMultiIndex(idx, r) {
  const out = new Array(r);
  for (let k = r - 1; k >= 0; k--) {
    out[k] = idx % 3;
    idx = Math.floor(idx / 3);
  }
  return out;
}
function applyVisualizerToggle() {
  const sw = document.getElementById("visualizerSwitch");
  const panel = document.getElementById("visualizerPanel");
  if (!sw || !panel) return;

  const apply = () => {
    const show = !!sw.checked;
    panel.style.display = show ? "block" : "none";

    // If we just opened it, ensure Three.js canvas resizes correctly.
    if (show) {
      requestAnimationFrame(() => {
        initViewer3D(); // safe: it resizes if already initialized :contentReference[oaicite:3]{index=3}

        const c = document.getElementById("viewer3d");
        if (c && renderer3d && camera3d) {
          const ww = c.clientWidth;
          const hh = c.clientHeight;
          if (ww && hh) {
            camera3d.aspect = ww / hh;
            camera3d.updateProjectionMatrix();
            renderer3d.setSize(ww, hh);
          }
        }

        const wSel = document.getElementById("wyckoffSelect");
        if (!wSel || !currentGroup || !wSel.value || wSel.value === "whole") {
          clearTensorSlidersUI();
        } else {
          updateOrbitViewer();
        }
      });
    }
  };

  sw.addEventListener("change", apply);
  apply(); // apply once on page load
}

// document.addEventListener("DOMContentLoaded", applyVisualizerToggle);
if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", applyVisualizerToggle);
} else {
  applyVisualizerToggle(); // DOM already loaded
}
function syncRankIsoWithMaxRank(maxRank) {
  const select = document.getElementById("rankIsoSelect");
  if (!select) return;

  const options = Array.from(select.options);
  let highestEnabled = null;

  options.forEach(opt => {
    const val = parseInt(opt.value, 10);
    if (Number.isNaN(val)) return;

    if (val <= maxRank) {
      opt.disabled = false;
      opt.hidden = false;
      highestEnabled = opt;
    } else {
      opt.disabled = true;
      opt.hidden = true;
    }
  });

  const currentVal = parseInt(select.value, 10);
  if (!Number.isNaN(currentVal) && currentVal <= maxRank) return;

  if (highestEnabled) {
    select.value = highestEnabled.value;
    select.dispatchEvent(new Event("change", { bubbles: true }));
  }
}


// === Tensor-Visualizer style: Tensor -> (allowed_harmonics, weights) -> sliders ===
// Performance tuning (lower values = faster, higher = smoother)
let VIS_NTH = 30;
let VIS_NPH = 60;
const PROJ_NTH = 36;
const PROJ_NPH = 72;

let _isoRenderHandle = null;
function scheduleIsoRender(rank, basis, spin = ACTIVE_SPIN) {
  if (_isoRenderHandle) cancelAnimationFrame(_isoRenderHandle);
  _isoRenderHandle = requestAnimationFrame(() => {
    _isoRenderHandle = null;
    renderTensorIsosurfacesFromHarmonics(rank, basis, spin);
  });
}
// Uses projectTensorToRealYlmLocal() which must implement the SAME convention as tensor_visualizer_copy.html:
// - integration Nth=60, Nph=120
// - consider all l <= rank
// - keep only non-zero harmonics (with eps ~ 1e-4)
// ==================== Tensor_visualizer-style projection (global) ====================
// Returns { harmonics: [[l,m],...], weights: [w...] }
// Expects: getRealYlm(l,m,theta,phi) exists
// Uses: _evalTensorOnDir(tensor, rank, nx, ny, nz) if available; otherwise uses fallback contraction.
function exportViewerPNG(filename = "scene.png", scale = 2) {
  if (!renderer3d || !scene3d || !camera3d) {
    alert("3D viewer not initialized yet.");
    return;
  }

  // Render one clean frame right now
  if (controls3d) controls3d.update();

  const container = document.getElementById("viewer3d");
  const w = container?.clientWidth || renderer3d.domElement.width;
  const h = container?.clientHeight || renderer3d.domElement.height;

  // Save renderer state
  const prevPixelRatio = renderer3d.getPixelRatio();
  const prevSize = renderer3d.getSize(new THREE.Vector2());

  // Higher-res export
  renderer3d.setPixelRatio((window.devicePixelRatio || 1) * scale);
  renderer3d.setSize(w, h, false);
  renderer3d.render(scene3d, camera3d);

  // Convert canvas -> PNG and download
  const dataURL = renderer3d.domElement.toDataURL("image/png");

  // Restore renderer state
  renderer3d.setPixelRatio(prevPixelRatio);
  renderer3d.setSize(prevSize.x, prevSize.y, false);
  renderer3d.render(scene3d, camera3d);

  const a = document.createElement("a");
  a.href = dataURL;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
}



// Fallback spatial contraction if your app.js doesn't already define _evalTensorOnDir
if (typeof window._evalTensorOnDir !== "function") {
  window._evalTensorOnDir = function(tensorFlat, rank, nx, ny, nz) {
    const n = [nx, ny, nz];
    let s = 0;
    for (let idx = 0; idx < tensorFlat.length; idx++) {
      const T = tensorFlat[idx];
      if (!T) continue;
      const mi = _idxToMultiIndex(idx, rank);
      let prod = 1;
      for (let k = 0; k < rank; k++) prod *= n[mi[k]];
      s += T * prod;
    }
    return s;
  };
}

window.projectTensorToRealYlmLocal = function projectTensorToRealYlmLocal(rank, tensorFlat, opts = {}) {
  const eps = opts.eps ?? 1e-4;
  const Nth = opts.Nth ?? PROJ_NTH;
  const Nph = opts.Nph ?? PROJ_NPH;
  const lMin = Number.isFinite(opts.lMin) ? opts.lMin : 0;
  const lMax = Number.isFinite(opts.lMax) ? opts.lMax : rank;

  // Scan harmonics l = lMin..lMax
  const harmonics = [];
  for (let l = lMin; l <= lMax; l++) {
    for (let m = -l; m <= l; m++) harmonics.push([l, m]);
  }

  const weights = solveHarmonicWeightsByScalarProduct(
    rank,
    tensorFlat,
    harmonics,
    { Nth, Nph }
  );

  // Filter zeros to reduce sliders (same eps logic)
  const outH = [];
  const outW = [];
  for (let i = 0; i < harmonics.length; i++) {
    if (Math.abs(weights[i]) > eps) {
      outH.push(harmonics[i]);
      outW.push(weights[i]);
    }
  }

  return { harmonics: outH, weights: outW };
};

// Compute harmonic weights using a scalar-product-aware fit (normal equations)
// to account for non-orthogonality introduced by finite sampling.
function solveHarmonicWeightsByScalarProduct(rank, tensorFlat, harmonics, opts = {}) {
  const Nth = opts.Nth ?? PROJ_NTH;
  const Nph = opts.Nph ?? PROJ_NPH;
  const nH = harmonics.length;

  if (nH === 0) return [];

  const rhs = new Array(nH).fill(0);
  const gram = Array.from({ length: nH }, () => new Array(nH).fill(0));

  const RES_THETA = Nth;
  const RES_PHI   = Nph;
  const dPhi   = (2 * Math.PI) / RES_PHI;
  const dTheta = Math.PI / RES_THETA;

  for (let t = 0; t < RES_THETA; t++) {
    const theta = (t + 0.5) * dTheta;
    const sinT = Math.sin(theta);
    const cosT = Math.cos(theta);
    const dOmega = sinT * dTheta * dPhi;

    for (let p = 0; p < RES_PHI; p++) {
      const phi = p * dPhi;

      const nx = sinT * Math.cos(phi);
      const ny = sinT * Math.sin(phi);
      const nz = cosT;

      const field = window._evalTensorOnDir(tensorFlat, rank, nx, ny, nz);

      const y = new Array(nH);
      for (let h = 0; h < nH; h++) {
        const [l, m] = harmonics[h];
        y[h] = getRealYlm(l, m, theta, phi);
      }

      for (let i = 0; i < nH; i++) {
        rhs[i] += field * y[i] * dOmega;
        for (let j = 0; j < nH; j++) {
          gram[i][j] += y[i] * y[j] * dOmega;
        }
      }
    }
  }

  for (let i = 0; i < nH; i++) gram[i][i] += 1e-8;

  try {
    const inv = invertMatrixN(gram);
    return matVec(inv, rhs);
  } catch (e) {
    console.warn("Harmonic fit failed, falling back to raw inner products.", e);
    return rhs;
  }
}

// Cache: per (fullRank, spin) store allowed harmonics for the whole tensor space
const ALLOWED_HARM_CACHE = {}; // key -> [[l,m],...]
const WIGNER_ROT_CACHE = new Map();

function _hmKey(l, m) { return `${l},${m}`; }

function buildFullHarmonics(maxL) {
  const out = [];
  for (let l = 0; l <= maxL; l++) {
    for (let m = -l; m <= l; m++) out.push([l, m]);
  }
  return out;
}

function filterHarmonicsByWeight(harmonics, weights, eps = 1e-6) {
  const outH = [];
  const outW = [];
  for (let i = 0; i < harmonics.length; i++) {
    if (Math.abs(weights[i]) > eps) {
      outH.push(harmonics[i]);
      outW.push(weights[i]);
    }
  }
  return { harmonics: outH, weights: outW };
}

function realYlmVector(l, theta, phi) {
  const vec = [];
  for (let m = -l; m <= l; m++) vec.push(getRealYlm(l, m, theta, phi));
  return vec;
}

function rotateDirection(R, v) {
  return [
    R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
    R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
    R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2]
  ];
}

function cartRotationMatrixFromFrac(R) {
  if (!R) return [[1,0,0],[0,1,0],[0,0,1]];
  const { a, b, c } = CELL_BASIS || { a:[1,0,0], b:[0,1,0], c:[0,0,1] };
  const B = new THREE.Matrix3().set(
    a[0], b[0], c[0],
    a[1], b[1], c[1],
    a[2], b[2], c[2]
  );
  const Binv = B.clone();
  if (Math.abs(Binv.determinant()) < 1e-12) return [[1,0,0],[0,1,0],[0,0,1]];
  Binv.invert();

  const Rf = new THREE.Matrix3().set(
    R[0][0], R[0][1], R[0][2],
    R[1][0], R[1][1], R[1][2],
    R[2][0], R[2][1], R[2][2]
  );
  const Rc = B.clone().multiply(Rf).multiply(Binv);
  const e = Rc.elements;
  return [
    [e[0], e[1], e[2]],
    [e[3], e[4], e[5]],
    [e[6], e[7], e[8]]
  ];
}

function matMul3x3(A, B) {
  const out = [[0,0,0],[0,0,0],[0,0,0]];
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      out[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j];
    }
  }
  return out;
}

function getRealYlmRotationMatrix(l, Rcart) {
  const key = `${l}|${Rcart.flat().map(v => v.toFixed(6)).join(',')}`;
  if (WIGNER_ROT_CACHE.has(key)) return WIGNER_ROT_CACHE.get(key);

  const dim = 2 * l + 1;
  const Nth = 16;
  const Nph = 32;
  const dTheta = Math.PI / Nth;
  const dPhi = (2 * Math.PI) / Nph;

  const Yorig = Array.from({ length: dim }, () => []);
  const Yrot = Array.from({ length: dim }, () => []);

  for (let t = 0; t < Nth; t++) {
    const theta = (t + 0.5) * dTheta;
    const sinT = Math.sin(theta);
    const cosT = Math.cos(theta);
    for (let p = 0; p < Nph; p++) {
      const phi = p * dPhi;
      const x = sinT * Math.cos(phi);
      const y = sinT * Math.sin(phi);
      const z = cosT;
      const v = [x, y, z];

      const vRot = rotateDirection(Rcart, v);
      const rxy = Math.sqrt(vRot[0] * vRot[0] + vRot[1] * vRot[1]);
      const thetaR = Math.acos(Math.max(-1, Math.min(1, vRot[2])));
      let phiR = Math.atan2(vRot[1], vRot[0]);
      if (phiR < 0) phiR += 2 * Math.PI;

      const yOrig = realYlmVector(l, theta, phi);
      const yRot = realYlmVector(l, thetaR, phiR);
      for (let i = 0; i < dim; i++) {
        Yorig[i].push(yOrig[i]);
        Yrot[i].push(yRot[i]);
      }
    }
  }

  const YorigT = matT(Yorig);
  const denom = matMul(Yorig, YorigT);
  for (let i = 0; i < dim; i++) denom[i][i] += 1e-8;
  const inv = invertMatrixN(denom);
  const numer = matMul(Yrot, YorigT);
  const D = matMul(numer, inv);

  WIGNER_ROT_CACHE.set(key, D);
  return D;
}

function rotateHarmonicWeightsByWigner(harmonics, weights, Rcart) {
  const byL = new Map();
  for (let i = 0; i < harmonics.length; i++) {
    const l = harmonics[i][0];
    if (!byL.has(l)) byL.set(l, []);
    byL.get(l).push({ idx: i, m: harmonics[i][1] });
  }

  const out = new Array(weights.length).fill(0);
  for (const [l, entries] of byL.entries()) {
    entries.sort((a, b) => a.m - b.m);
    const dim = 2 * l + 1;
    if (entries.length !== dim) {
      // if incomplete, keep weights unchanged for this l
      for (const e of entries) out[e.idx] = weights[e.idx];
      continue;
    }
    const w = entries.map(e => weights[e.idx]);
    const D = getRealYlmRotationMatrix(l, Rcart);
    const wRot = matVec(D, w);
    for (let i = 0; i < entries.length; i++) {
      out[entries[i].idx] = wRot[i];
    }
  }
  return out;
}

// Project ONLY on a provided harmonic list (no filtering, stable slider set)
function projectTensorToSpecificHarmonics(rank, tensorFlat, harmonics, opts = {}) {
  const Nth = opts.Nth ?? PROJ_NTH;
  const Nph = opts.Nph ?? PROJ_NPH;

  return solveHarmonicWeightsByScalarProduct(
    rank,
    tensorFlat,
    harmonics,
    { Nth, Nph }
  );
}

// Discover allowed harmonics from the WHOLE basis space (union)
// This is the Option A behavior = tensor_visualizer variable discovery.
function discoverAllowedHarmonicsFromBasis(fullRank, basis, spinIndex, eps = 1e-4) {
  const key = `${fullRank}_s${spinIndex}`;

  if (ALLOWED_HARM_CACHE[key]) return ALLOWED_HARM_CACHE[key];

  const active = new Set();

  for (const B of basis) {
    const { spatialRank, spatialTensor } = getSpatialTensorForProjection(
      fullRank,
      spinIndex,
      B,
      (CELL_BASIS && CELL_BASIS.a) ? [
        [CELL_BASIS.a[0], CELL_BASIS.b[0], CELL_BASIS.c[0]],
        [CELL_BASIS.a[1], CELL_BASIS.b[1], CELL_BASIS.c[1]],
        [CELL_BASIS.a[2], CELL_BASIS.b[2], CELL_BASIS.c[2]]
      ] : [[1,0,0],[0,1,0],[0,0,1]]
    );

    // We can use your existing projector with filtering eps to speed up.
    // It returns only non-zero harmonics for this basis vector.
    const proj = projectTensorToRealYlmLocal(
      spatialRank,
      new Float64Array(spatialTensor),
      { eps, Nth: PROJ_NTH, Nph: PROJ_NPH, lMin: spatialRank, lMax: spatialRank }
    );

    for (const [l, m] of proj.harmonics) active.add(_hmKey(l, m));
  }

  // Convert set -> sorted list by (l then m)
  const out = Array.from(active)
    .map(s => s.split(',').map(Number))
    .sort((a,b) => (a[0]-b[0]) || (a[1]-b[1]));

  ALLOWED_HARM_CACHE[key] = out;
  return out;
}
// ===== tensor_visualizer-style VARIABLE persistence (here: basis coefficients) =====
// rank(fullRank) -> coeffs array (length = basis.length)
const GLOBAL_BASIS_COEFFS = {};   // like GLOBAL_VARIABLE_VALUES in tensor_visualizer
// Registry is declared at top-of-file to keep labels stable across atoms

let TV_STATE = {}; // rank -> { harmonics: [[l,m],...], weights: [...] }

// rank = FULL rank coming from backend (includes spin index as first dimension)
// fullRank includes spin index as first dimension
function buildSiteYlmSliders(fullRank, basis) {
  ensureSpinButtons();

  const box = document.getElementById('tensor-sliders');
  const forbidden = document.getElementById('tensor-forbidden');
  if (!box) return;

  // keep spin row, clear rest
  const spinRow = document.getElementById('site-spin-row');
  box.innerHTML = '';
  if (spinRow) box.appendChild(spinRow);

  if (!basis || basis.length === 0) {
    if (forbidden) {
      forbidden.style.display = 'block';
      forbidden.innerText = 'No basis for this rank.';
    }
    return;
  }

  // === Discover allowed harmonics from WHOLE basis space ===
  const allowed = discoverAllowedHarmonicsFromBasis(fullRank, basis, ACTIVE_SPIN, 1e-4);

  if (allowed.length === 0) {
    if (forbidden) {
      forbidden.style.display = 'block';
      forbidden.innerText = 'No active harmonics found for this rank/spin.';
    }
    return;
  }
  if (forbidden) forbidden.style.display = 'none';

  // --- helper: compare harmonic lists ---
  const sameHarmonics = (A, B) => {
    if (!A || !B || A.length !== B.length) return false;
    for (let i = 0; i < A.length; i++) {
      if (A[i][0] !== B[i][0] || A[i][1] !== B[i][1]) return false;
    }
    return true;
  };

  // --- key: memory is per (rank, spin) ---
  const key = `${fullRank}_s${ACTIVE_SPIN}`;

  // Ensure GLOBAL_BASIS_COEFFS exists for this rank (initialized to zeros)
  if (!GLOBAL_BASIS_COEFFS[fullRank] || GLOBAL_BASIS_COEFFS[fullRank].length !== basis.length) {
    GLOBAL_BASIS_COEFFS[fullRank] = new Array(basis.length).fill(0);
  }

  if (!TV_STATE[key] || !sameHarmonics(TV_STATE[key].harmonics, allowed)) {
    TV_STATE[key] = {
      harmonics: allowed,
      weights: new Array(allowed.length).fill(0)
    };
    TV_STATE[key].sliderEls = [];
  }

  // Populate harmonic weights from global basis coefficients
  if (GLOBAL_BASIS_COEFFS[fullRank] && GLOBAL_BASIS_COEFFS[fullRank].length === basis.length) {
    try {
      const P = buildBasisToHarmonicMatrix(fullRank, basis, ACTIVE_SPIN, allowed);
      const W = matVec(P, GLOBAL_BASIS_COEFFS[fullRank]);
      TV_STATE[key].weights = W.map(x => (Math.abs(x) < 1e-6 ? 0 : x));
    } catch (e) {
      console.warn("Spin inherit failed (basis->harmonics)", e);
    }
  }

  const orbitalLabel = (lVal, mVal) => {
    const labels = {
      0: { 0: "s" },
      1: { "-1": "p<sub>y</sub>", 0: "p<sub>z</sub>", 1: "p<sub>x</sub>" },
      2: {
        "-2": "d<sub>xy</sub>",
        "-1": "d<sub>yz</sub>",
        0: "d<sub>z<sup>2</sup></sub>",
        1: "d<sub>xz</sub>",
        2: "d<sub>x<sup>2</sup>-y<sup>2</sup></sub>"
      },
      3: {
        "-3": "f<sub>y(3x<sup>2</sup>-y<sup>2</sup>)</sub>",
        "-2": "f<sub>xyz</sub>",
        "-1": "f<sub>yz<sup>2</sup></sub>",
        0: "f<sub>z<sup>3</sup></sub>",
        1: "f<sub>xz<sup>2</sup></sub>",
        2: "f<sub>z(x<sup>2</sup>-y<sup>2</sup>)</sub>",
        3: "f<sub>x(x<sup>2</sup>-3y<sup>2</sup>)</sub>"
      },
      4: {
        "-4": "g<sub>xy(x<sup>2</sup>-y<sup>2</sup>)</sub>",
        "-3": "g<sub>y(3x<sup>2</sup>-y<sup>2</sup>)z</sub>",
        "-2": "g<sub>xyz<sup>2</sup></sub>",
        "-1": "g<sub>yz<sup>3</sup></sub>",
        0: "g<sub>z<sup>4</sup></sub>",
        1: "g<sub>xz<sup>3</sup></sub>",
        2: "g<sub>z<sup>2</sup>(x<sup>2</sup>-y<sup>2</sup>)</sub>",
        3: "g<sub>x(x<sup>2</sup>-3y<sup>2</sup>)z</sub>",
        4: "g<sub>x<sup>4</sup>+y<sup>4</sup></sub>"
      }
    };

    const label = labels[lVal]?.[mVal] ?? labels[lVal]?.[String(mVal)];
    return label || `l=${lVal}, m=${mVal}`;
  };

  // === Build sliders for ALL allowed harmonics ===
  let currentL = null;
  for (let i = 0; i < allowed.length; i++) {
    const [l, m] = allowed[i];

    if (l !== currentL) {
      const h = document.createElement('div');
      h.innerHTML = `<b>L = ${l}</b>`;
      h.style.marginTop = "10px";
      h.style.borderBottom = "1px solid #eee";
      box.appendChild(h);
      currentL = l;
    }

    const row = document.createElement('div');
    row.style.display = 'flex';
    row.style.alignItems = 'center';
    row.style.gap = '8px';
    row.style.margin = '4px 0';

    const label = document.createElement('div');
    label.style.width = '160px';
    label.style.fontFamily = 'monospace';
    label.innerHTML = `Y(${l},${m}) ${orbitalLabel(l, m)}`;

    const slider = document.createElement('input');
    slider.type = 'range';
    slider.min = -2;
    slider.max =  2;
    slider.step = 0.05;

    slider.value = TV_STATE[key].weights[i];

    const value = document.createElement('div');
    value.style.width = '55px';
    value.style.fontFamily = 'monospace';
    value.innerText = Number(slider.value).toFixed(2);

    slider.oninput = () => {
      const v = parseFloat(slider.value);
      TV_STATE[key].weights[i] = v;
      value.innerText = v.toFixed(2);

      try {
        const P_active = buildBasisToHarmonicMatrix(fullRank, basis, ACTIVE_SPIN, TV_STATE[key].harmonics);
        const c = leastSquares(P_active, TV_STATE[key].weights, 1e-8);
        GLOBAL_BASIS_COEFFS[fullRank] = c;
        scheduleIsoRender(fullRank, basis, ACTIVE_SPIN);
      } catch (e) {
        console.warn("Global basis sync failed:", e);
      }
    };

    row.appendChild(label);
    row.appendChild(slider);
    row.appendChild(value);
    box.appendChild(row);
    TV_STATE[key].sliderEls.push({ slider, value, l, m });
  }
}



function tensorSignedValueFromWeights(harmonics, weights, theta, phi) {
  let s = 0;
  for (let i = 0; i < harmonics.length; i++) {
    const w = weights[i];
    if (Math.abs(w) < 1e-6) continue;
    const [l, m] = harmonics[i];
    s += w * getRealYlm(l, m, theta, phi);
  }
  return s; // IMPORTANT: keep sign
}


function tensorRadiusFromWeights(harmonics, weights, theta, phi) {
  let s = 0;
  for (let i = 0; i < harmonics.length; i++) {
    const w = weights[i];
    if (Math.abs(w) < 1e-6) continue;
    const [l, m] = harmonics[i];
    s += w * getRealYlm(l, m, theta, phi);
  }
  return Math.abs(s);
}

function buildBlobGeometryFromWeights(harmonics, weights, scale = 0.14) {
  const Nth = VIS_NTH;
  const Nph = VIS_NPH;

  const positions = [];
  const indices = [];
  const svals = [];   // signed values per vertex
  const vid = (i, j) => i * Nph + j;

  // --- build vertices + store signed values ---
  for (let i = 0; i < Nth; i++) {
    const theta = Math.PI * i / (Nth - 1);
    for (let j = 0; j < Nph; j++) {
      const phi = 2 * Math.PI * j / (Nph - 1);

      const s = tensorSignedValueFromWeights(harmonics, weights, theta, phi);
      const r = Math.abs(s);
      
      svals.push(s);

      const x = Math.sin(theta) * Math.cos(phi);
      const y = Math.sin(theta) * Math.sin(phi);
      const z = Math.cos(theta);

      positions.push(scale * r * x, scale * r * y, scale * r * z);
    }
  }

  // --- faces ---
  for (let i = 0; i < Nth - 1; i++) {
    for (let j = 0; j < Nph - 1; j++) {
      const a = vid(i, j);
      const b = vid(i + 1, j);
      const c = vid(i + 1, j + 1);
      const d = vid(i, j + 1);
      indices.push(a, b, d);
      indices.push(b, c, d);
    }
  }

  // --- colors: pure sign mapping (negative vs positive) ---
  const colors = [];
  const mode = document.getElementById('multipoleTypeSelect')?.value || 'magnetic';
  const posColor = (mode === 'electric') ? [0.98, 0.62, 0.12] : [1.0, 0.0, 0.0];
  const negColor = (mode === 'electric') ? [0.12, 0.62, 0.68] : [0.0, 0.0, 1.0];
  for (let k = 0; k < svals.length; k++) {
    const isNeg = svals[k] < 0;
    const [r, g, b] = isNeg ? negColor : posColor;
    colors.push(r, g, b);
  }

  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
  geo.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3)); // <--- KEY
  geo.setIndex(indices);
  geo.computeVertexNormals();
  return geo;
}
function resetWyckoffSelectionAndClearViewer() {
  const sel = document.getElementById("wyckoffSelect");
  if (sel) {
    sel.value = "";
    sel.selectedIndex = 0;
  }

  wyckoffRandomVars = null;

  // clear the 3D scene
  clearViewerOrbit();
  clearIsoGroup();

  // clear cached orbit used by tensor rendering
  _lastOrbitForTensor = null;

  // ✅ BUG 5 FIX: clear stale sliders from last position
  clearTensorSlidersUI();
}

function renderTensorIsosurfacesFromHarmonics(rank, basis, spin = ACTIVE_SPIN) {
  if (!scene3d) return;

  const orbit = _lastOrbitForTensor;
  if (!orbit) return;

  if (!basis || basis.length === 0) {
    clearIsoGroup();
    return;
  }

  if (!GLOBAL_BASIS_COEFFS[rank] || GLOBAL_BASIS_COEFFS[rank].length !== basis.length) {
    GLOBAL_BASIS_COEFFS[rank] = new Array(basis.length).fill(0);
  }

  const coeffs = GLOBAL_BASIS_COEFFS[rank];
  const fullRank = rank;
  const spatialRank = Math.max(0, fullRank - 1);
  const dim = Math.pow(3, fullRank);
  const tmpBuf = new Float64Array(dim);

  if (isoGroup3d) scene3d.remove(isoGroup3d);
  isoGroup3d = new THREE.Group();

  const orbitExpanded = expandOrbitForViewer(orbit);
  const { a, b, c: cc } = CELL_BASIS || { a:[1,0,0], b:[0,1,0], c:[0,0,1] };
  
  // DEBUG BASIS
  console.log("DEBUG: CELL_BASIS a=", a, "b=", b);

  const M = [ [a[0], b[0], cc[0]], [a[1], b[1], cc[1]], [a[2], b[2], cc[2]] ];


  for (const atom of orbitExpanded) {
    const [u, v, w] = atom.coord;
    const [x, y, z] = fracToCart(u, v, w);

    const R_frac = atom.op?.R;
    const tr = (atom.op?.tr !== undefined) ? atom.op.tr : 1;

    const localFrac = new Float64Array(dim);
    for (let bIdx = 0; bIdx < basis.length; bIdx++) {
      const coeff = coeffs[bIdx];
      if (Math.abs(coeff) < 1e-12) continue;
      const baseVec = basis[bIdx].slice();
      const transformed = R_frac ? applyOpToTensor(baseVec, fullRank, R_frac, tmpBuf) : baseVec;
      for (let i = 0; i < dim; i++) localFrac[i] += coeff * transformed[i];
    }

  const mMode = (typeof multipoleMode !== "undefined") ? multipoleMode : 
          (document.getElementById('multipoleTypeSelect')?.value || 'electric');

    if (mMode === "magnetic") {
      const detR = R_frac ? det3x3(R_frac) : 1;
      const signFactor = tr * detR;
      if (Math.abs(signFactor - 1) > 1e-12) {
        for (let i = 0; i < localFrac.length; i++) localFrac[i] *= signFactor;
      }
    }

    const { spatialRank: projRank, spatialTensor } = getSpatialTensorForProjection(
      fullRank,
      spin,
      localFrac,
      M
    );

    const proj = projectTensorToRealYlmLocal(projRank, spatialTensor, {
      eps: 1e-4,
      Nth: PROJ_NTH,
      Nph: PROJ_NPH,
      lMin: projRank,
      lMax: projRank
    });

    // DEBUG: Monitor Y00 weights
    if (proj.harmonics) {
       const idx00 = proj.harmonics.findIndex(h => h[0] === 0 && h[1] === 0);
       if (idx00 !== -1) {
         const label = (atom.originalIndex !== undefined) ? `Atom ${atom.originalIndex + 1}` : 'Atom ?';
         console.log(`[Y00 Monitor] ${label} at (${u.toFixed(3)}, ${v.toFixed(3)}, ${w.toFixed(3)}) : ${proj.weights[idx00].toPrecision(5)}`);
       }
    }

    if (!proj.harmonics || proj.harmonics.length === 0) continue;

    const geo = buildBlobGeometryFromWeights(proj.harmonics, proj.weights, 0.14);
    const mesh = new THREE.Mesh(
      geo,
      new THREE.MeshStandardMaterial({
        vertexColors: true,
        transparent: false,
        opacity: 1.0,
        side: THREE.DoubleSide
      })
    );
    mesh.renderOrder = 1;
    mesh.position.set(x, y, z);
    isoGroup3d.add(mesh);
  }

  scene3d.add(isoGroup3d);
}

let _lastOrbitForTensor = null;
// Returns the currently selected tensor rank from the UI
function getSelectedRank() {
  const el = document.getElementById("rankIsoSelect");
  if (!el) return null;
  const val = parseInt(el.value, 10);
  return Number.isNaN(val) ? null : val;
}

// Return the basis array for a given rank from the latest backend results.
// lastResults entries look like: {rank: number, basis: [...], group, wyckoff}
function getBasisForRank(rank) {
  if (!Array.isArray(lastResults)) return [];
  const r = lastResults.find(x => x && x.rank === rank);
  return (r && Array.isArray(r.basis)) ? r.basis : [];
}
function renderTensorIsosurfacesFromSiteTensor(orbitMaybe) {
  if (!scene3d) return;

  if (orbitMaybe) _lastOrbitForTensor = orbitMaybe;
  const orbit = _lastOrbitForTensor;
  if (!orbit) return;

  // Clear old group
  if (isoGroup3d) scene3d.remove(isoGroup3d);
  isoGroup3d = new THREE.Group();

  const rank = getSelectedRank();
  if (!rank) return;

    const basis = getBasisForRank(rank);
    if (!basis || basis.length === 0) {
      clearIsoGroup();
      // clear any stale slider state for this rank/spin
      const key = `${rank}_s${ACTIVE_SPIN}`;
      if (TV_STATE && TV_STATE[key]) delete TV_STATE[key];
      buildSiteYlmSliders(rank, basis);
      return;
    }

    buildSiteYlmSliders(rank, basis);

  renderTensorIsosurfacesFromHarmonics(rank, basis, ACTIVE_SPIN);



}

// Bootstrap
        loadGroupsFromBackend().catch(e => {
            console.error(e);
            document.getElementById('status').innerText = 'Error loading groups: ' + e.message;
        });
document.getElementById('wyckoffSelect').onchange = async () => {
  renderWyckoffOccupancyLine();

  // Full reset on Wyckoff change
  resetAllForWyckoffChange();

  wyckoffRandomVars = null;
  await triggerCompute();
  updateOrbitViewer();
};


const rankIsoSelect = document.getElementById("rankIsoSelect");
if (rankIsoSelect) {
  rankIsoSelect.addEventListener("change", () => {
    resetTensorSliderMemory();   //  Bug 6 rule
    clearTensorSlidersUI();      // optional
    updateOrbitViewer();
  });
}

// buildYlmSlidersUI() is called after each compute() so sliders match symmetry
document.getElementById("exportSceneBtn").onclick = () => {
  exportViewerPNG("tensor_scene.png", 2); // 2 = nice quality
};





