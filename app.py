"""
PTM Site Scanner - 100% Pure Python Standard Library Web App (Group 15).
Zero Streamlit, zero Flask, and zero external HTML/CSS/JS files.
Run in VS Code: python app.py  ->  open http://127.0.0.1:8000
"""
from email.parser import BytesParser
from email.policy import default
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import os
from urllib.parse import parse_qs, quote, urlparse
import xml.etree.ElementTree as ET

from ptm_data import MOD_INFO, MOTIFS
from scanner import scan_sequences
from sequence_reader import read_sequences_from_text

SAMPLE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "samples")

BUILTIN_SAMPLES = {
    "Demo FASTA (p53 Kinase & SUMO Motifs)": (
        ">sp|P04637|P53_HUMAN Cellular tumor antigen p53\n"
        "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP\n"
        "DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK\n"
        "SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE\n"
        "RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS\n"
        "SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP\n"
        "PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG\n"
        "GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD\n"
    ),
    "Demo PDB (SEQRES Chain A)": (
        "SEQRES   1 A   25  MET SER PRO ARG ARG LYS SER GLY ASN ALA THR CYS\n"
        "SEQRES   2 A   25  VAL LYS GLU ASP PRO GLY ARG GLY PHE TYR SER THR\n"
        "SEQRES   3 A   25  PRO\n"
    ),
}


def get_available_samples() -> dict:
    """Loads any files in ./samples plus built-in FASTA/PDB examples."""
    samples = dict(BUILTIN_SAMPLES)
    if os.path.isdir(SAMPLE_DIR):
        for fname in sorted(os.listdir(SAMPLE_DIR)):
            fpath = os.path.join(SAMPLE_DIR, fname)
            if os.path.isfile(fpath):
                with open(fpath, encoding="utf-8", errors="replace") as f:
                    samples[fname] = f.read()
    return samples


def parse_multipart_form(content_type: str, body: bytes):
    """Parses file uploads and form fields using Python's built-in email parser."""
    header = f"Content-Type: {content_type}\r\nMIME-Version: 1.0\r\n\r\n".encode("utf-8")
    msg = BytesParser(policy=default).parsebytes(header + body)
    fields = {}
    files = {}
    if msg.is_multipart():
        for part in msg.iter_parts():
            name = part.get_param("name", header="content-disposition")
            filename = part.get_param("filename", header="content-disposition")
            payload = part.get_payload(decode=True) or b""
            if filename:
                files[name] = (filename, payload)
            elif name:
                fields[name] = payload.decode("utf-8", errors="replace")
    return fields, files


class PTMWebPage:
    """Constructs the PTM Scanner web interface purely via Python ElementTree."""

    def __init__(self, title: str):
        self.root = ET.Element("html")
        head = ET.SubElement(self.root, "head")
        ET.SubElement(head, "meta", charset="utf-8")
        t = ET.SubElement(head, "title")
        t.text = title

        body = ET.SubElement(self.root, "body", bgcolor="#f1f5f9", text="#0f172a")

        # Header banner
        banner = ET.SubElement(
            body, "table", width="100%", bgcolor="#0f172a", cellpadding="16", cellspacing="0"
        )
        b_tr = ET.SubElement(banner, "tr")
        b_td1 = ET.SubElement(b_tr, "td")
        f_title = ET.SubElement(b_td1, "font", face="Arial, sans-serif", size="5", color="#ffffff")
        b_bold = ET.SubElement(f_title, "b")
        b_bold.text = "PTM Site Scanner — Group 15"

        ET.SubElement(b_td1, "br")
        f_sub = ET.SubElement(b_td1, "font", face="Arial, sans-serif", size="2", color="#94a3b8")
        f_sub.text = "Post-Translational Modification & Motif Predictor (100% Pure Python)"

        b_td2 = ET.SubElement(b_tr, "td", align="right")
        home_link = ET.SubElement(b_td2, "a", href="/")
        f_home = ET.SubElement(home_link, "font", face="Arial, sans-serif", size="3", color="#38bdf8")
        f_home_b = ET.SubElement(f_home, "b")
        f_home_b.text = "New Scan / Home"

        # Main centered layout container
        self.container = ET.SubElement(
            body, "table", width="960", align="center", cellpadding="10", cellspacing="0"
        )

    def add_card(self, heading: str, subtitle: str = None, bg: str = "#ffffff", border: str = "#cbd5e1"):
        row = ET.SubElement(self.container, "tr")
        cell = ET.SubElement(row, "td")
        card = ET.SubElement(
            cell,
            "table",
            width="100%",
            bgcolor=bg,
            border="1",
            bordercolor=border,
            cellpadding="16",
            cellspacing="0",
        )
        tr = ET.SubElement(card, "tr")
        td = ET.SubElement(tr, "td")

        h_font = ET.SubElement(td, "font", face="Arial, sans-serif", size="4", color="#0f172a")
        h_bold = ET.SubElement(h_font, "b")
        h_bold.text = heading

        if subtitle:
            p = ET.SubElement(td, "p")
            s_font = ET.SubElement(p, "font", face="Arial, sans-serif", size="2", color="#475569")
            s_font.text = subtitle

        return td

    def add_table(self, parent_td, headers, rows):
        tbl = ET.SubElement(
            parent_td,
            "table",
            width="100%",
            bgcolor="#ffffff",
            border="1",
            bordercolor="#cbd5e1",
            cellpadding="8",
            cellspacing="0",
        )
        h_tr = ET.SubElement(tbl, "tr", bgcolor="#e2e8f0")
        for h in headers:
            th = ET.SubElement(h_tr, "th", align="left")
            f = ET.SubElement(th, "font", face="Arial, sans-serif", size="2", color="#0f172a")
            f.text = str(h)

        for r in rows:
            tr = ET.SubElement(tbl, "tr")
            for val in r:
                td = ET.SubElement(tr, "td", valign="top")
                if isinstance(val, ET.Element):
                    td.append(val)
                else:
                    f = ET.SubElement(td, "font", face="Arial, sans-serif", size="2", color="#1e293b")
                    f.text = str(val)

    def add_scanner_form(self, prefill_text: str = "", lenient_checked: bool = False):
        td = self.add_card(
            "Scan Protein Sequence(s)",
            "Paste a raw amino acid sequence, FASTA record(s), or PDB SEQRES lines, or upload a file (.fasta, .fa, .faa, .pdb, .ent).",
        )

        # Sample loader links
        samples = get_available_samples()
        sample_p = ET.SubElement(td, "p")
        lbl = ET.SubElement(sample_p, "font", face="Arial, sans-serif", size="2", color="#334155")
        lbl_b = ET.SubElement(lbl, "b")
        lbl_b.text = "Load Example: "

        for idx, sname in enumerate(samples.keys()):
            if idx > 0:
                sep = ET.SubElement(sample_p, "font", face="Arial, sans-serif", size="2", color="#94a3b8")
                sep.text = "  |  "
            a = ET.SubElement(sample_p, "a", href=f"/?sample={quote(sname)}")
            af = ET.SubElement(a, "font", face="Arial, sans-serif", size="2", color="#0284c7")
            af.text = sname

        # Upload & sequence form
        form = ET.SubElement(
            td, "form", method="POST", action="/scan", enctype="multipart/form-data"
        )

        p_file = ET.SubElement(form, "p")
        f_lbl = ET.SubElement(p_file, "font", face="Arial, sans-serif", size="2", color="#0f172a")
        f_lbl_b = ET.SubElement(f_lbl, "b")
        f_lbl_b.text = "1. Upload Sequence File (optional): "
        ET.SubElement(
            p_file,
            "input",
            type="file",
            name="file",
            accept=".pdb,.ent,.fasta,.fa,.faa,.txt",
        )

        p_txt = ET.SubElement(form, "p")
        t_lbl = ET.SubElement(p_txt, "font", face="Arial, sans-serif", size="2", color="#0f172a")
        t_lbl_b = ET.SubElement(t_lbl, "b")
        t_lbl_b.text = "2. Or Paste Sequence / FASTA / PDB Text:"
        ET.SubElement(p_txt, "br")
        ta = ET.SubElement(p_txt, "textarea", name="text", rows="9", cols="95")
        ta.text = prefill_text or ""

        p_opt = ET.SubElement(form, "p")
        cb_attrs = {"type": "checkbox", "name": "lenient", "value": "true"}
        if lenient_checked:
            cb_attrs["checked"] = "checked"
        ET.SubElement(p_opt, "input", **cb_attrs)
        opt_lbl = ET.SubElement(p_opt, "font", face="Arial, sans-serif", size="2", color="#334155")
        opt_lbl.text = " Lenient mode (replace unknown/non-standard characters with X instead of raising an error)"

        p_sub = ET.SubElement(form, "p")
        ET.SubElement(p_sub, "input", type="submit", value="Run PTM Scan")

    def add_reference_guide(self):
        td = self.add_card(
            "Supported Modification Rules & Consensus Motifs",
            "Reference lookup rules and regular-expression motifs evaluated by scanner.py:",
        )
        motif_rows = [
            [m["modification"], m["name"], m["pattern"], MOD_INFO.get(m["modification"], "")]
            for m in MOTIFS
        ]
        self.add_table(
            td,
            ["Modification", "Motif Rule", "Regex Pattern", "Biological Description"],
            motif_rows,
        )

    def render(self) -> bytes:
        html_str = ET.tostring(self.root, encoding="unicode", method="html")
        return ("<!DOCTYPE html>\n" + html_str).encode("utf-8")


class PTMRequestHandler(BaseHTTPRequestHandler):
    def address_string(self):
        # Bypass slow reverse-DNS lookups on incoming requests
        return self.client_address[0]

    def send_html(self, page: PTMWebPage, status: int = 200):
        data = page.render()
        self.send_response(status)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        parsed = urlparse(self.path)
        if parsed.path == "/favicon.ico":
            self.send_response(204)
            self.end_headers()
            return

        if parsed.path == "/scan":
            self.send_response(303)
            self.send_header("Location", "/")
            self.end_headers()
            return

        if parsed.path != "/":
            page = PTMWebPage("404 Not Found")
            page.add_card("404 — Page Not Found", "Click 'New Scan / Home' above to return.")
            self.send_html(page, status=404)
            return

        params = parse_qs(parsed.query)
        sample_key = params.get("sample", [None])[0]
        prefill = ""
        if sample_key:
            prefill = get_available_samples().get(sample_key, "")

        page = PTMWebPage("PTM Site Scanner | Group 15")
        page.add_scanner_form(prefill_text=prefill)
        page.add_reference_guide()
        self.send_html(page)

    def do_POST(self):
        parsed = urlparse(self.path)
        if parsed.path != "/scan":
            self.send_error(404, "Not Found")
            return

        content_length = int(self.headers.get("Content-Length", 0))
        content_type = self.headers.get("Content-Type", "")
        raw_body = self.rfile.read(content_length)

        fields, files = parse_multipart_form(content_type, raw_body)
        lenient = fields.get("lenient") == "true"
        text_input = fields.get("text", "")
        uploaded = files.get("file")

        if uploaded and uploaded[0]:
            filename = uploaded[0]
            content = uploaded[1].decode("utf-8", errors="replace")
        elif text_input.strip():
            filename, content = None, text_input
        else:
            page = PTMWebPage("Input Error | PTM Site Scanner")
            page.add_card(
                "No Sequence Provided",
                "Please upload a FASTA/PDB file or paste a sequence before running the scan.",
                bg="#fef2f2",
                border="#fca5a5",
            )
            page.add_scanner_form(prefill_text=text_input, lenient_checked=lenient)
            self.send_html(page, status=400)
            return

        try:
            sequences = read_sequences_from_text(content, filename)
        except ValueError as err:
            page = PTMWebPage("File Format Error | PTM Site Scanner")
            page.add_card("File Parsing Error", str(err), bg="#fef2f2", border="#fca5a5")
            page.add_scanner_form(prefill_text=text_input, lenient_checked=lenient)
            self.send_html(page, status=400)
            return

        if not sequences:
            page = PTMWebPage("Empty Input | PTM Site Scanner")
            page.add_card(
                "No Sequences Found",
                "PDB files must contain SEQRES lines, and FASTA files must start with a '>' header line.",
                bg="#fef2f2",
                border="#fca5a5",
            )
            page.add_scanner_form(prefill_text=text_input, lenient_checked=lenient)
            self.send_html(page, status=400)
            return

        report = scan_sequences(sequences, lenient=lenient)
        page = PTMWebPage("Scan Results | PTM Site Scanner")

        for chain in report["chains"]:
            cname = chain["name"]
            if "error" in chain:
                page.add_card(
                    f"Chain / Record: {cname} — Error",
                    chain["error"],
                    bg="#fef2f2",
                    border="#fca5a5",
                )
                continue

            td = page.add_card(
                f"Chain / Record: {cname} (Length: {chain['length']} aa)",
                f"Total candidate PTM sites detected: {len(chain['sites'])}",
            )

            for warn in chain["warnings"]:
                wp = ET.SubElement(td, "p")
                wf = ET.SubElement(wp, "font", face="Arial, sans-serif", size="2", color="#b45309")
                wb = ET.SubElement(wf, "b")
                wb.text = f"Warning: {warn}"

            seq_p = ET.SubElement(td, "p")
            seq_lbl = ET.SubElement(seq_p, "font", face="Arial, sans-serif", size="2", color="#334155")
            seq_lbl_b = ET.SubElement(seq_lbl, "b")
            seq_lbl_b.text = "Cleaned Sequence: "
            code_el = ET.SubElement(seq_p, "code")
            code_el.text = chain["sequence"]

            if chain["counts"]:
                summary_rows = []
                for mod_name, stat in chain["counts"].items():
                    summary_rows.append([
                        mod_name,
                        stat["total"],
                        stat["motif"],
                        MOD_INFO.get(mod_name, ""),
                    ])
                page.add_table(
                    td,
                    [
                        "Modification Type",
                        "Total Candidate Sites",
                        "Motif-Matched Sites",
                        "Biological Description",
                    ],
                    summary_rows,
                )
                ET.SubElement(td, "br")

            if chain["sites"]:
                site_rows = []
                for site in chain["sites"]:
                    for mod in site["modifications"]:
                        if mod["motifs"]:
                            motif_desc = "; ".join(
                                f"{m['name']} [{m['match']} at {m['start']}-{m['end']}]"
                                for m in mod["motifs"]
                            )
                            evidence_label = "Motif Match (High Confidence)"
                        else:
                            motif_desc = "Residue-level rule only"
                            evidence_label = "Residue Lookup"

                        site_rows.append([
                            site["position"],
                            f"{site['residue']} ({site['residue_name']})",
                            mod["name"],
                            evidence_label,
                            motif_desc,
                        ])

                page.add_table(
                    td,
                    ["Position", "Residue", "Predicted PTM", "Evidence", "Matched Motif & Window"],
                    site_rows,
                )
            else:
                no_p = ET.SubElement(td, "p")
                no_f = ET.SubElement(no_p, "font", face="Arial, sans-serif", size="2", color="#475569")
                no_f.text = "No modifiable residues found in this sequence."

        page.add_scanner_form(prefill_text=content if not uploaded else "", lenient_checked=lenient)
        self.send_html(page)


if __name__ == "__main__":
    host = "0.0.0.0"
    port = int(os.environ.get("PORT", 8000))
    server = ThreadingHTTPServer((host, port), PTMRequestHandler)
    server.daemon_threads = True
    print(f"PTM Site Scanner live! Click to open: http://127.0.0.1:{port}")
    print("Press Ctrl+C to stop the server.")
    server.serve_forever()