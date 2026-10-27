// Gene lookup for the NPM1 AML project page.
// Reads data/genes.json (built by scripts/12_build_site.R). Each gene row is
// [symbol, ensembl_id, biotype, expr, de, note, tcga]:
//   expr = 3 groups x [lower whisker, q1, median, q3, upper whisker]
//   de   = Beat AML [log2FC, lfcSE, padj, rank, significant (0/1)] or null
//   tcga = TCGA-LAML [log2FC, lfcSE, padj, rank, significant (0/1),
//          symbol if renamed since Ensembl 75 or null] or null (not tested)
(function () {
  "use strict";

  var COLORS = ["#D55E00", "#0072B2", "#8c8c8c"];
  var data = null;          // parsed genes.json
  var bySymbol = {};        // UPPERCASE symbol -> [row index, ...]
  var byId = {};            // Ensembl ID -> row index
  var symbolsSorted = [];   // unique uppercase symbols, sorted, for prefix search

  var q = document.getElementById("q");
  var list = document.getElementById("suggest");
  var statusEl = document.getElementById("status");
  var result = document.getElementById("result");
  var active = -1;          // highlighted suggestion

  function esc(s) {
    return String(s).replace(/[&<>"']/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
    });
  }

  function fmtP(p) {
    if (p === null || p === undefined) return "not available";
    if (p === 0) return "&lt; 1e-300";
    if (p < 0.001) {
      var e = p.toExponential(1).split("e");
      return e[0] + " &times; 10<sup>" + Number(e[1]) + "</sup>";
    }
    return p.toPrecision(2);
  }

  function fmtNum(x, d) { return (x > 0 ? "+" : "") + x.toFixed(d); }

  // ---- Data loading -------------------------------------------------------
  fetch("data/genes.json")
    .then(function (r) {
      if (!r.ok) throw new Error("HTTP " + r.status);
      return r.json();
    })
    .then(function (json) {
      data = json;
      data.genes.forEach(function (g, i) {
        var k = g[0].toUpperCase();
        (bySymbol[k] = bySymbol[k] || []).push(i);
        byId[g[1].toUpperCase()] = i;
      });
      symbolsSorted = Object.keys(bySymbol).sort();
      q.disabled = false;
      statusEl.textContent = data.genes.length.toLocaleString() + " genes loaded.";
      var fromHash = decodeURIComponent((location.hash.match(/gene=([^&]+)/) || [])[1] || "");
      showByQuery(fromHash || "HOXA9", !fromHash);
    })
    .catch(function (err) {
      statusEl.textContent = "Couldn't load the gene data (" + err.message + "). " +
        "If you opened this file directly from disk, serve the docs/ folder instead, " +
        "for example with: python3 -m http.server --directory docs";
    });

  // ---- Search -----------------------------------------------------------
  function matches(text) {
    var t = text.trim().toUpperCase();
    if (!t) return [];
    var out = [];
    if (/^ENSG\d+/.test(t)) {
      if (byId[t] !== undefined) out.push(byId[t]);
      return out;
    }
    // Exact symbol first, then prefix matches in alphabetical order
    (bySymbol[t] || []).forEach(function (i) { out.push(i); });
    var lo = 0, hi = symbolsSorted.length;
    while (lo < hi) { var mid = (lo + hi) >> 1; if (symbolsSorted[mid] < t) lo = mid + 1; else hi = mid; }
    for (var j = lo; j < symbolsSorted.length && out.length < 8; j++) {
      var s = symbolsSorted[j];
      if (s.indexOf(t) !== 0) break;
      if (s === t) continue;
      bySymbol[s].forEach(function (i) { if (out.length < 8) out.push(i); });
    }
    return out;
  }

  function renderSuggest() {
    var m = data ? matches(q.value) : [];
    active = -1;
    if (!m.length) { list.hidden = true; q.setAttribute("aria-expanded", "false"); return; }
    list.innerHTML = m.map(function (i, k) {
      var g = data.genes[i];
      return '<li role="option" id="opt' + k + '" data-i="' + i + '">' +
        "<span>" + esc(g[0]) + "</span><small>" + esc(g[1]) + "</small></li>";
    }).join("");
    list.hidden = false;
    q.setAttribute("aria-expanded", "true");
  }

  function setActive(k) {
    var items = list.querySelectorAll("li");
    if (!items.length) return;
    active = (k + items.length) % items.length;
    items.forEach(function (li, n) { li.setAttribute("aria-selected", n === active ? "true" : "false"); });
    q.setAttribute("aria-activedescendant", "opt" + active);
    items[active].scrollIntoView({ block: "nearest" });
  }

  function closeSuggest() { list.hidden = true; q.setAttribute("aria-expanded", "false"); }

  q.addEventListener("input", renderSuggest);
  q.addEventListener("keydown", function (e) {
    if (e.key === "ArrowDown") { e.preventDefault(); if (list.hidden) renderSuggest(); setActive(active + 1); }
    else if (e.key === "ArrowUp") { e.preventDefault(); setActive(active - 1); }
    else if (e.key === "Escape") { closeSuggest(); }
    else if (e.key === "Enter") {
      e.preventDefault();
      var items = list.querySelectorAll("li");
      if (!list.hidden && active >= 0 && items[active]) show(Number(items[active].dataset.i));
      else showByQuery(q.value);
    }
  });
  list.addEventListener("mousedown", function (e) {
    var li = e.target.closest("li");
    if (li) { e.preventDefault(); show(Number(li.dataset.i)); }
  });
  q.addEventListener("blur", function () { setTimeout(closeSuggest, 100); });

  document.querySelectorAll(".chips button").forEach(function (b) {
    b.addEventListener("click", function () { showByQuery(b.dataset.gene); });
  });

  window.addEventListener("hashchange", function () {
    var g = decodeURIComponent((location.hash.match(/gene=([^&]+)/) || [])[1] || "");
    if (g && data) showByQuery(g, true);
  });

  function showByQuery(text, keepHash) {
    if (!data) return;
    var m = matches(text);
    var t = text.trim().toUpperCase();
    var exact = m.filter(function (i) { return data.genes[i][0].toUpperCase() === t || data.genes[i][1] === t; });
    if (exact.length) show(exact[0], keepHash);
    else if (m.length) { renderSuggest(); statusEl.textContent = "No exact match for “" + text.trim() + "”. Pick one of the suggestions."; }
    else {
      closeSuggest();
      statusEl.textContent = "“" + text.trim() + "” isn't in the data. It may not pass the low-count filter, " +
        "or it may have an older symbol in Ensembl 75; try its Ensembl ID.";
    }
  }

  // ---- Result -----------------------------------------------------------
  function show(i, keepHash) {
    var g = data.genes[i];
    var sym = g[0], id = g[1], biotype = g[2], expr = g[3], de = g[4], note = g[5], tcga = g[6];
    closeSuggest();
    q.value = sym;
    if (!keepHash) history.replaceState(null, "", "#gene=" + encodeURIComponent(sym));

    var dupes = bySymbol[sym.toUpperCase()].length;
    statusEl.textContent = dupes > 1
      ? sym + " matches " + dupes + " Ensembl IDs; showing " + id + ". Search by ID to see the others."
      : "";

    document.getElementById("g-symbol").textContent = sym;
    document.getElementById("g-meta").innerHTML =
      '<a href="https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=' + encodeURIComponent(id) +
      '" rel="noopener">' + esc(id) + "</a>" + (biotype ? " &middot; " + esc(biotype.replace(/_/g, " ")) : "");
    document.getElementById("g-note").innerHTML = note ? '<div class="note">' + esc(note) + "</div>" : "";

    // Unhide before drawing: Plotly measures the box width when it draws, and a
    // hidden box measures 0, which makes it fall back to a fixed 700 px chart.
    result.hidden = false;
    drawPlot(expr);

    var med = [0, 1, 2].map(function (k) { return expr[k * 5 + 2]; });
    document.getElementById("g-expr-text").textContent =
      "Median VST: NPM1-mutated " + med[0].toFixed(1) + ", KMT2A-rearranged " + med[1].toFixed(1) +
      ", other AML " + med[2].toFixed(1) + ".";

    document.getElementById("g-de").innerHTML = deTable(de, tcga);
  }

  // ---- DESeq2 results, Beat AML next to TCGA-LAML --------------------------
  // One column per cohort. A cohort without a result gets one merged cell
  // with the reason instead of five empty ones.
  var DE_ROWS = ["Call", "log2FC &plusmn; SE", "Fold change", "padj (BH)", "Rank"];

  function deCells(r, nGenes) {
    var lfc = r[0], se = r[1], padj = r[2], rank = r[3], sig = r[4] === 1;
    var call = sig ? (lfc > 0 ? '<span class="badge up">Higher</span>' : '<span class="badge down">Lower</span>')
                   : '<span class="badge ns">Not changed</span>';
    var fold = Math.pow(2, Math.abs(lfc));
    var foldTxt = Math.abs(lfc) < 0.05 ? "about the same"
      : fold.toFixed(fold < 10 ? 2 : 1) + "x " + (lfc > 0 ? "higher" : "lower");
    return [call, fmtNum(lfc, 2) + " &plusmn; " + se.toFixed(2), foldTxt, fmtP(padj),
            rank.toLocaleString() + " of " + nGenes.toLocaleString()];
  }

  function deTable(de, tcga) {
    var cols = [
      de ? deCells(de, data.meta.de_n_genes)
         : "Not in the model: passes the low-count filter across all 300 samples, but not across the 276 main-cohort samples.",
      tcga ? deCells(tcga, data.meta.tcga_n_genes)
           : "Not tested in TCGA: no matching Ensembl ID in GENCODE v36, or too few reads for TCGA's low-count filter."
    ];
    var body = DE_ROWS.map(function (label, r) {
      var cells = cols.map(function (c) {
        if (typeof c !== "string") return "<td>" + c[r] + "</td>";
        return r === 0 ? '<td class="missing" rowspan="' + DE_ROWS.length + '">' + esc(c) + "</td>" : "";
      }).join("");
      return '<tr><th scope="row">' + label + "</th>" + cells + "</tr>";
    }).join("");
    var renamed = tcga && tcga[5]
      ? '<p class="small">In TCGA (GENCODE v36) this gene is called <b>' + esc(tcga[5]) + "</b>; matched on the Ensembl ID.</p>"
      : "";
    return '<table class="de-table"><thead><tr><td></td>' +
      '<th scope="col">Beat AML 2.0<span>n = 276 (100 vs 176)</span></th>' +
      '<th scope="col">TCGA-LAML<span>n = 111 (31 vs 80)</span></th>' +
      "</tr></thead><tbody>" + body + "</tbody></table>" + renamed;
  }

  function cssVar(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }

  function drawPlot(expr) {
    var el = document.getElementById("plot");
    if (typeof Plotly === "undefined") {
      el.style.height = "auto";  // no empty chart-sized gap above the numbers
      el.innerHTML = '<p class="msg">The plotting library (Plotly, from cdn.plot.ly) didn\'t load. The numbers below still apply.</p>';
      return;
    }
    var n = data.meta.n;
    var short = ["NPM1-mutated", "KMT2A-r", "Other AML"];  // short labels fit a phone screen
    var traces = data.meta.groups.map(function (grp, k) {
      var s = expr.slice(k * 5, k * 5 + 5);
      var label = short[k] + "<br>n = " + n[k];
      return {
        type: "box", name: label, x: [label],
        lowerfence: [s[0]], q1: [s[1]], median: [s[2]], q3: [s[3]], upperfence: [s[4]],
        marker: { color: COLORS[k] }, line: { color: COLORS[k], width: 2 },
        fillcolor: COLORS[k] + "55", boxpoints: false, hoverinfo: "y+name"
      };
    });
    var text = cssVar("--text"), line = cssVar("--line"), surface = cssVar("--surface");
    var layout = {
      yaxis: { title: { text: "Expression (VST, ~log2)" }, zeroline: false, gridcolor: line, color: text },
      xaxis: { color: text, fixedrange: true, tickangle: 0 },
      showlegend: false,
      margin: { l: 60, r: 10, t: 16, b: 56 },
      paper_bgcolor: surface, plot_bgcolor: surface,
      font: { family: "-apple-system, BlinkMacSystemFont, Segoe UI, Roboto, Arial, sans-serif", color: text }
    };
    Plotly.react(el, traces, layout, { responsive: true, displayModeBar: false });
  }
})();
