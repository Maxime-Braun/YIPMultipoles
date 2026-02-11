// MCIF and Wyckoff parsing utilities extracted from main.js

export function parseMcifValue(raw) {
  if (!raw) return null;
  const cleaned = String(raw)
    .replace(/["']/g, "")
    .replace(/\([^)]*\)/g, "")
    .trim();
  const val = parseFloat(cleaned);
  return Number.isFinite(val) ? val : null;
}

export function parseMcifString(raw) {
  if (!raw) return null;
  return String(raw).replace(/[\"']/g, "").trim();
}

export function tokenizeMcifRow(line) {
  if (!line) return [];
  const tokens = line.match(/'(?:[^']*)'|"(?:[^"]*)"|\S+/g);
  return tokens ? tokens.map(t => t.trim()) : [];
}

export function extractMcifRaw(text, tag) {
  const re = new RegExp(`${tag}\\s+([^\\r\\n]+)`, "i");
  const match = text.match(re);
  return match ? match[1].trim() : null;
}

export function extractMcifTag(text, tag) {
  const raw = extractMcifRaw(text, tag);
  return raw ? parseMcifValue(raw) : null;
}

export function parseMcifMagneticGroup(text) {
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

export function parseMcifAtomFract(text) {
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

export function parseMcifAtomCoords(text) {
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
      const lowerHeaders = headers.map(h => h.toLowerCase().trim());
      const labelIdx = lowerHeaders.findIndex(h =>
        h === "_atom_site_label" ||
        h === "_atom_site_type_symbol" ||
        h === "_atom_site_symbol" ||
        h.includes("_atom_site_label") ||
        h.includes("_atom_site_type_symbol")
      );
      const xIdx = lowerHeaders.findIndex(h => h === "_atom_site_fract_x");
      const yIdx = lowerHeaders.findIndex(h => h === "_atom_site_fract_y");
      const zIdx = lowerHeaders.findIndex(h => h === "_atom_site_fract_z");
      // Add wyckoff extraction
      const wyckoffIdx = lowerHeaders.findIndex(h =>
        h === "_atom_site_wyckoff_symbol" ||
        h === "_atom_site_wyckoff_letter" ||
        h === "_atom_site_wyckoff"
      );
      const multIdx = lowerHeaders.findIndex(h => h === "_atom_site_symmetry_multiplicity");
      if (labelIdx >= 0 && xIdx >= 0 && yIdx >= 0 && zIdx >= 0) {
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
            const x = parseMcifValue(parts[xIdx]);
            const y = parseMcifValue(parts[yIdx]);
            const z = parseMcifValue(parts[zIdx]);
            let wyckoff = "";
            if (wyckoffIdx >= 0) {
              const wyckoffSym = parseMcifString(parts[wyckoffIdx]);
              let mult = "";
              if (multIdx >= 0) {
                const multRaw = parseMcifString(parts[multIdx]);
                mult = multRaw ? String(multRaw).replace(/\D/g, "") : "";
              }
              wyckoff = wyckoffSym ? `${mult}${wyckoffSym}` : "";
            }
            if (label && [x, y, z].every(v => v !== null)) {
              sites.push({ label, coord: [x, y, z], wyckoff });
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

export function parseMcifMagneticLabels(text) {
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
      const labelIdx = lowerHeaders.findIndex(h => h === "_atom_site_moment.label");
      if (labelIdx >= 0) {
        const labels = [];
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
            if (label) labels.push(label);
          }
          i += 1;
        }
        return labels;
      }
    }
    i += 1;
  }
  return [];
}

export function parseMcifAtomSites(text) {
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
      const lowerHeaders = headers.map(h => h.toLowerCase().trim());
      // Build a map from header name to index for flexible lookup
      const headerMap = {};
      lowerHeaders.forEach((h, idx) => { headerMap[h] = idx; });
      // Try to find the best match for each field
      function findHeaderIdx(possibles) {
        for (const p of possibles) {
          const idx = lowerHeaders.findIndex(h => h === p || h.includes(p));
          if (idx >= 0) return idx;
        }
        return -1;
      }
      const labelIdx = findHeaderIdx([
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_symbol"
      ]);
      const wyckoffIdx = findHeaderIdx([
        "_atom_site_wyckoff_symbol",
        "_atom_site_wyckoff_letter",
        "_atom_site_wyckoff"
      ]);
      const multIdx = findHeaderIdx([
        "_atom_site_symmetry_multiplicity"
      ]);
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
          // Use headerMap to get the right value, even if columns are missing or extra
          const label = labelIdx < parts.length ? parseMcifString(parts[labelIdx]) : null;
          const wyckoffSym = wyckoffIdx < parts.length ? parseMcifString(parts[wyckoffIdx]) : null;
          const multRaw = (multIdx >= 0 && multIdx < parts.length) ? parseMcifString(parts[multIdx]) : null;
          const mult = multRaw ? String(multRaw).replace(/\D/g, "") : "";
          const wyckoff = wyckoffSym ? `${mult}${wyckoffSym}` : null;
          if (label && wyckoff) {
            sites.push({ label, wyckoff });
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

export function renderMcifAtomHeader(atomSites) {
  const header = document.getElementById("mcifAtomHeader");
  if (!header) return;
  if (!atomSites || atomSites.length === 0) {
    header.innerHTML = "";
    header.style.display = "none";
    return;
  }
  const formatLabel = (label) => {
    if (!label) return "";
    const match = String(label).match(/^(.*?)(\d+)$/);
    if (!match) return String(label);
    const base = match[1];
    const digits = match[2];
    return `${base}<sub>${digits}</sub>`;
  };
  const magneticLabels = Array.isArray(window.MCIF_MAGNETIC_LABELS)
    ? new Set(window.MCIF_MAGNETIC_LABELS.map(l => String(l).toLowerCase()))
    : new Set();
  const items = atomSites
    .map(site => {
      const isMagnetic = site && site.label && magneticLabels.has(String(site.label).toLowerCase());
      const labelHtml = formatLabel(site.label);
      const labelSpan = isMagnetic
        ? `<span style="color:#c62828; font-weight:700;">${labelHtml}</span>`
        : labelHtml;
      // Always use the wyckoff site from the site object if present
      const wyckoff = site && site.wyckoff ? site.wyckoff : "";
      return `${labelSpan}: ${wyckoff}`;
    })
    .join(", ");
  header.innerHTML = `Atoms (Wyckoff): ${items}`;
  header.style.display = "block";
}

export function parseWyckoffComponent(s) {
  const res = { x: 0, y: 0, z: 0, c: 0 };
  if (!s) return res;
  let expr = s.replace(/\s+/g, '').toLowerCase();
  let terms = expr.replace(/-/g, '+-').split('+').filter(t => t !== '');
  terms.forEach(t => {
    let coeff = 1;
    if (t.startsWith('-')) {
      coeff = -1;
      t = t.substring(1);
    }
    if (t.startsWith('+')) t = t.substring(1);
    if (t.includes('x')) {
      let val = t.replace('x', '');
      if (val === '') val = '1';
      res.x += coeff * parseFloat(val);
    } else if (t.includes('y')) {
      let val = t.replace('y', '');
      if (val === '') val = '1';
      res.y += coeff * parseFloat(val);
    } else if (t.includes('z')) {
      let val = t.replace('z', '');
      if (val === '') val = '1';
      res.z += coeff * parseFloat(val);
    } else if (t) {
      if (t.includes('/')) {
        const [n, d] = t.split('/');
        res.c += coeff * (parseFloat(n) / parseFloat(d));
      } else {
        res.c += coeff * parseFloat(t);
      }
    }
  });
  return res;
}

export function evalWyckoffComponent(comp, vars) {
  return comp.x * vars.x + comp.y * vars.y + comp.z * vars.z + comp.c;
}

export function normalizeFrac(val) {
  let v = val - Math.floor(val + 1e-8);
  if (Math.abs(v - 1.0) < 1e-6) v = 0.0;
  return v;
}

export function fracClose(a, b, tol = 1e-3) {
  const diff = Math.abs(a - b);
  return diff < tol || Math.abs(diff - 1) < tol;
}

export function solveWyckoffVars(components, targetCoord) {
  const A = components.map(c => [c.x, c.y, c.z]);
  const b = components.map((c, i) => normalizeFrac(targetCoord[i] - c.c));
  const det = window.det3x3 ? window.det3x3(A) : 0;
  if (Math.abs(det) > 1e-6 && window.invertMatrix3x3) {
    const shifts = [-1, 0, 1];
    for (let dx of shifts) for (let dy of shifts) for (let dz of shifts) {
      const bShift = [b[0] + dx, b[1] + dy, b[2] + dz];
      const inv = window.invertMatrix3x3(A);
      const vars = {
        x: inv[0][0] * bShift[0] + inv[0][1] * bShift[1] + inv[0][2] * bShift[2],
        y: inv[1][0] * bShift[0] + inv[1][1] * bShift[1] + inv[1][2] * bShift[2],
        z: inv[2][0] * bShift[0] + inv[2][1] * bShift[1] + inv[2][2] * bShift[2]
      };
      if ([vars.x, vars.y, vars.z].every(v => isFinite(v))) {
        return vars;
      }
    }
  }
  const vars = { x: 0, y: 0, z: 0 };
  const solved = { x: false, y: false, z: false };
  components.forEach((c, idx) => {
    const coeffs = [c.x, c.y, c.z];
    const nonZero = coeffs.filter(v => Math.abs(v) > 1e-6).length;
    if (nonZero === 1) {
      const val = normalizeFrac(targetCoord[idx] - c.c);
      if (Math.abs(c.x) > 1e-6) {
        vars.x = val / c.x;
        solved.x = true;
      } else if (Math.abs(c.y) > 1e-6) {
        vars.y = val / c.y;
        solved.y = true;
      } else if (Math.abs(c.z) > 1e-6) {
        vars.z = val / c.z;
        solved.z = true;
      }
    }
  });
  if (solved.x || solved.y || solved.z) return vars;
  return null;
}

export function inferWyckoffLabelFromCoord(group, coord) {
  if (!group || !Array.isArray(group.wyckoff)) return null;
  const matches = [];
  for (const w of group.wyckoff) {
    if (!w || !w.coord) continue;
    const parts = w.coord.split(',');
    if (parts.length !== 3) continue;
    const comps = parts.map(parseWyckoffComponent);
    const vars = solveWyckoffVars(comps, coord);
    let isMatch = false;
    if (!vars) {
      const fixed = comps.every((c, idx) =>
        Math.abs(c.x) < 1e-6 && Math.abs(c.y) < 1e-6 && Math.abs(c.z) < 1e-6 &&
        fracClose(normalizeFrac(c.c), normalizeFrac(coord[idx]))
      );
      isMatch = fixed;
    } else {
      const evalCoord = comps.map(c => normalizeFrac(evalWyckoffComponent(c, vars)));
      const matchesCoord = evalCoord.every((v, idx) => fracClose(v, normalizeFrac(coord[idx])));
      isMatch = matchesCoord;
    }
    if (isMatch) {
      const label = w.label || '';
      const multMatch = label.match(/^(\d+)/);
      const multiplicity = multMatch ? parseInt(multMatch[1], 10) : Number.POSITIVE_INFINITY;
      const varCount = comps.reduce((count, c) => {
        const varsHere = [c.x, c.y, c.z].filter(v => Math.abs(v) > 1e-6).length;
        return count + (varsHere > 0 ? 1 : 0);
      }, 0);
      matches.push({ label, multiplicity, varCount });
    }
  }
  if (!matches.length) return null;
  matches.sort((a, b) => {
    if (a.multiplicity !== b.multiplicity) return a.multiplicity - b.multiplicity;
    return a.varCount - b.varCount;
  });
  return matches[0].label || null;
}

export function getMcifOccupiedMagneticSites() {
  if (!window.MCIF_ATOM_SITES || window.MCIF_ATOM_SITES.length === 0) return null;
  // If magnetic labels are provided, restrict to those atom labels
  const magneticLabels = Array.isArray(window.MCIF_MAGNETIC_LABELS) ? window.MCIF_MAGNETIC_LABELS.map(l => String(l).toLowerCase()) : null;
  let labels = [];
  if (magneticLabels && magneticLabels.length) {
    const setMag = new Set(magneticLabels);
    labels = window.MCIF_ATOM_SITES
      .filter(site => site && site.label && setMag.has(String(site.label).toLowerCase()))
      .map(site => site.wyckoff)
      .filter(Boolean)
      .map(v => String(v).trim())
      .filter(Boolean);
  } else {
    labels = window.MCIF_ATOM_SITES
      .map(site => site.wyckoff)
      .filter(Boolean)
      .map(val => String(val).trim())
      .filter(Boolean);
  }
  if (labels.length === 0) return null;
  return Array.from(new Set(labels));
}

export function getMcifMagneticSitesFromCoords() {
  if (!window.currentGroup || !Array.isArray(window.MCIF_ATOM_COORDS) || window.MCIF_ATOM_COORDS.length === 0) {
    return null;
  }
  const magneticLabels = Array.isArray(window.MCIF_MAGNETIC_LABELS) ? window.MCIF_MAGNETIC_LABELS : [];
  if (magneticLabels.length === 0) return null;
  const labelSet = new Set(magneticLabels.map(l => String(l).toLowerCase()));
  const inferred = [];
  window.MCIF_ATOM_COORDS.forEach(site => {
    if (!site || !site.label || !Array.isArray(site.coord)) return;
    if (!labelSet.has(String(site.label).toLowerCase())) return;
    const wyck = inferWyckoffLabelFromCoord(window.currentGroup, site.coord);
    if (wyck) inferred.push(wyck);
  });
  return inferred.length ? Array.from(new Set(inferred)) : null;
}

export function buildMcifAtomSitesWithWyckoff(group) {
  // Always infer wyckoff for all atoms with coords, regardless of MCIF_ATOM_SITES content
  if (!group || !Array.isArray(window.MCIF_ATOM_COORDS) || window.MCIF_ATOM_COORDS.length === 0) return window.MCIF_ATOM_SITES || [];
  console.debug('[DEBUG buildMcifAtomSitesWithWyckoff] group:', group && group.name ? group.name : group, 'window.MCIF_ATOM_SITES:', window.MCIF_ATOM_SITES, 'window.MCIF_ATOM_COORDS:', window.MCIF_ATOM_COORDS);
  const inferred = window.MCIF_ATOM_COORDS.map(site => {
    const wyckoff = inferWyckoffLabelFromCoord(group, site.coord);
    if (!wyckoff) {
      console.debug('[DEBUG buildMcifAtomSitesWithWyckoff] could not infer wyckoff for site:', site);
    }
    return { label: site.label, wyckoff: wyckoff || "" };
  });
  console.debug('[DEBUG buildMcifAtomSitesWithWyckoff] inferred:', inferred);
  return inferred.length > 0 ? inferred : (window.MCIF_ATOM_SITES || []);
}

export function parseMcifCellParams(text) {
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
