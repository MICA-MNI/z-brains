#!/usr/bin/env python
import os
import pandas as pd

# Snakemake variables
input_files = snakemake.input
output_files = snakemake.output
params = snakemake.params

# Assign variables from snakemake
score_files = input_files.get('score_files', [])
if isinstance(score_files, str):
    score_files = [score_files]
demographics_csv = input_files['demographics'] if 'demographics' in input_files else params.get('demographics_csv', None)
output_dir = params.get('output_dir', None)
subject = params.get('subject', None)
session = params.get('session', None)
verbose = params.get('verbose', True)
report_file = output_files['report'] if 'report' in output_files else output_files[0]

# Lookup age/sex from demographics
if demographics_csv is None:
    raise ValueError("demographics_csv must be provided via snakemake.input or snakemake.params")
df = pd.read_csv(demographics_csv)
row = df[(df['ID'] == subject) & (df['SES'] == session)]
if row.empty:
    raise ValueError(f"No entry for subject {subject}, session {session} in demographics CSV.")
age = row.iloc[0]['AGE'] if 'AGE' in row.columns else None
sex = row.iloc[0]['SEX'] if 'SEX' in row.columns else None
os.makedirs(output_dir, exist_ok=True)

# --- Report generation logic (as before) ---
def report_header_template(sid, ses=None, age=None, sex=None, analysis=None):
    return f"""
    <h1>Clinical Report for {sid} {f'({ses})' if ses else ''}</h1>
    <p>Age: {age if age is not None else 'N/A'} | Sex: {sex if sex is not None else 'N/A'} | Analysis: {analysis if analysis else 'N/A'}</p>
    <hr>
    """

def convert_html_to_pdf(source_html, output_filename):
    from xhtml2pdf import pisa
    with open(output_filename, "w+b") as result_file:
        pisa_status = pisa.CreatePDF(source_html, dest=result_file)
    return pisa_status.err

# For demonstration, just list all score files in the report
report_html = report_header_template(subject, ses=session, age=age, sex=sex)
report_html += "<h2>Included Score Files:</h2><ul>"
for score_file in score_files:
    report_html += f"<li>{score_file}</li>"
report_html += "</ul>"
# Example: add a placeholder figure
fig_path = os.path.join(output_dir, f"{subject}_example.png")
import matplotlib.pyplot as plt
import numpy as np
plt.figure(figsize=(4,2))
plt.plot(np.linspace(-3,3,100), np.tanh(np.linspace(-3,3,100)), label="Example Curve")
plt.title("Example Plot")
plt.xlabel("x")
plt.ylabel("tanh(x)")
plt.legend()
plt.tight_layout()
plt.savefig(fig_path)
plt.close()
report_html += f'<img src="{fig_path}" width="400"><br>'
report_html += "<p>Report generation complete. (Add more content as needed.)</p>"
convert_html_to_pdf(report_html, report_file)
if verbose:
    print(f"Report saved to {report_file}") 