# LD_X11: Reference table of the Green et al. 2024 astrocyte states (Ast.1-Ast.10).
# Content is a hand-curated summary of the source publication, not derived from the
# data, so this script has no inputs - it only renders the table.
# Styling matches the LD_X01 cohort tables (three rules, italic header, 9 pt).
# Outputs: CSV + styled HTML (always) and Word/PNG (flextable).

library(tidyverse)

### paths -------------------------------------------------------------------
# Canonical location is RDS on the HPC; also runs locally against the mounted
# share (renders in seconds, no data read) - the first existing path wins.
base_candidates = c("/rds/general/user/lvd25/home/AST_scRNAseq_TREM2",   # HPC
                    "/Volumes/lvd25/home/AST_scRNAseq_TREM2")            # RDS mounted locally
base = base_candidates[dir.exists(base_candidates)][1]
if (is.na(base)) stop("Neither RDS path is reachable - is the share mounted?")
out_dir = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X11_"
message("Using base: ", base, "\nWriting outputs to: ", out_dir)

### table content -----------------------------------------------------------
states = tribble(
  ~State,    ~Annotation,                          ~Processes, ~`Marker genes`,
  "Ast.1",   "Homeostatic-like",                   "Not specified", "Not named",
  "Ast.2",   "Homeostatic-like",                   "Not specified", "Not named",
  "Ast.3",   "Enhanced mitophagy and translation",
  "Translation, mitochondrial function, mitophagy, chaperone-mediated autophagy, oxidative phosphorylation",
  "PINK1",
  "Ast.4",   "Reactive-like",
  "Extracellular matrix organisation, excitatory synaptic genes", "GFAP, DPP10",
  "Ast.5",   "Reactive-like",
  "Axonogenesis, wound healing", "GFAP, SERPINA3, OSMR",
  "Ast.6",   "Not annotated in the source",        "Not specified", "Not named",
  "Ast.7",   "Interferon-responding",              "Interferon response", "IFI6",
  "Ast.8",   "Stress response",
  "Chemical and heat stress, sterol metabolism", "Not named",
  "Ast.9",   "Stress response",
  "Heat and oxidative stress, tau binding, necroptosis", "DNAJB1, HSPH1",
  "Ast.10",  "Stress response, AD-elevated",
  "Oxidative stress, reactive oxygen species, metallothioneins, zinc ion homeostasis",
  "SLC38A2, SMTN"
)

cap = paste("Supplementary Table 4. Astrocyte states defined by Green et al. (2024).",
            "Annotations, associated processes and marker genes are as reported in the",
            "source publication; states without a published annotation or named markers",
            "are indicated as such.")

### save CSV ----------------------------------------------------------------
write_csv(states, file.path(out_dir, paste0(script_ind, "Green24_astrocyte_states.csv")))

### render: styled HTML (no extra packages, opens in Word) ------------------
esc = function(x) { x = gsub("&", "&amp;", x); x = gsub("<", "&lt;", x); gsub(">", "&gt;", x) }

th = paste0("<th class='", c("s", "a", "p", "g"), "'>", esc(names(states)), "</th>", collapse = "")
body = map_chr(seq_len(nrow(states)), function(i) {
  cls = c("s", "a", "p", "g")
  cells = map_chr(seq_along(states),
                  ~ paste0("<td class='", cls[.x], "'>", esc(states[[.x]][i]), "</td>"))
  paste0("<tr>", paste(cells, collapse = ""), "</tr>")
})

css = paste(
  "body{font-family:'Times New Roman',Georgia,serif;font-size:10pt;margin:24px;}",
  "table{border-collapse:collapse;margin:0 0 28px 0;max-width:960px;}",
  "th,td{padding:4px 12px 4px 4px;text-align:left;vertical-align:top;}",
  "thead th{border-top:1.5px solid #000;border-bottom:1.5px solid #000;font-style:italic;font-weight:bold;}",
  "tbody tr:last-child td{border-bottom:1.5px solid #000;}",
  "td.s{font-weight:bold;white-space:nowrap;}",          # State
  "td.g{font-style:italic;}",                            # Marker genes (gene symbols)
  "td.p{max-width:420px;}",                              # Processes (wraps)
  "p.cap{font-size:9pt;margin:0 0 6px 0;max-width:900px;}", sep = "\n")

html_doc = paste0("<!DOCTYPE html>\n<html><head><meta charset='utf-8'>\n<style>\n", css,
                  "\n</style>\n</head><body>\n<p class='cap'>", esc(cap), "</p>\n",
                  "<table>\n<thead><tr>", th, "</tr></thead>\n<tbody>\n",
                  paste(body, collapse = "\n"), "\n</tbody>\n</table>\n</body></html>\n")

html_path = file.path(out_dir, paste0(script_ind, "Green24_astrocyte_states.html"))
con = file(html_path, open = "w", encoding = "UTF-8")
writeLines(html_doc, con); close(con)
message("Wrote styled HTML (open in Word): ", html_path)

### render (flextable -> Word/PNG) ------------------------------------------
if (!requireNamespace("flextable", quietly = TRUE)) {
  stop("flextable is not available to this R session. Check .libPaths() / the ",
       "loaded R module, or install it. The styled HTML above is already written ",
       "and can be used in the meantime.")
} else {
  library(flextable)

  say = function(label, expr) {
    ok = tryCatch({ force(expr); TRUE },
                  error = function(e) { message("  FAILED ", label, ": ", conditionMessage(e)); FALSE })
    if (ok) message("  wrote ", label)
    invisible(ok)
  }

  ft = flextable(states) %>%
    bold(part = "header") %>% italic(part = "header") %>%
    bold(j = "State", part = "body") %>%
    italic(j = "Marker genes", part = "body") %>%   # gene symbols
    valign(valign = "top", part = "body") %>%
    border_remove() %>%
    hline_top(part = "header", border = fp_border_default(width = 1.5)) %>%
    hline_bottom(part = "header", border = fp_border_default(width = 1.5)) %>%
    hline_bottom(part = "body", border = fp_border_default(width = 1.5)) %>%
    fontsize(size = 9, part = "all") %>%
    set_caption(cap) %>%
    width(j = c("State", "Annotation", "Processes", "Marker genes"),
          width = c(0.7, 1.9, 3.6, 1.5)) %>%
    line_spacing(space = 1.1, part = "body")

  message("flextable ", as.character(packageVersion("flextable")), " - rendering:")
  say("Green24_astrocyte_states.docx",
      save_as_docx(ft, path = file.path(out_dir, paste0(script_ind, "Green24_astrocyte_states.docx"))))
  # PNG needs webshot2/chromote (headless Chrome); the docx is the thesis deliverable
  say("Green24_astrocyte_states.png",
      save_as_image(ft, path = file.path(out_dir, paste0(script_ind, "Green24_astrocyte_states.png")), res = 300))
}

message("Done. Outputs in: ", out_dir)
