from xhtml2pdf import pisa
import sys
import pandas as pd
import os
import argparse

# defined command line options
CLI=argparse.ArgumentParser()
CLI.add_argument("subject", nargs=1, type=str, default="")
CLI.add_argument("session", nargs=1, type=str, default="")
CLI.add_argument("out", nargs=1, type=str, default="")
CLI.add_argument("demo", nargs=1, type=str, default="")
CLI.add_argument("thr", nargs=1, type=float, default=1.96)
CLI.add_argument(
  "--featList_ctx",  
  nargs="*",
  type=str,
  default=['flair', 'qt1', 'adc', 'thickness'])
CLI.add_argument(
  "--featList_sctx",
  nargs="*",
  type=str,  
  default=[])
CLI.add_argument(
  "--featList_hipp",
  nargs="*",
  type=str,  
  default=[])

# parse the command line
args = CLI.parse_args()

# access CLI options
subject = args.subject[0]
session = args.session[0]
out = args.out[0]
demo = args.demo[0]
thr = args.thr[0]
featList_ctx = ', '.join(args.featList_ctx)
featList_sctx = ', '.join(args.featList_sctx)
featList_hipp = ', '.join(args.featList_hipp)

TBL = pd.read_excel(demo, engine='openpyxl', skiprows=[0,2]).dropna(how="all").dropna(axis=1, how="all")

# Get participant's sex
if TBL[TBL['ID'] == subject]['sex'].to_list()[0] == 'F':
    sex = 'Female'
elif TBL[TBL['ID'] == subject]['sex'].to_list()[0] == 'M':
    sex = 'Male'
else:
    sex = 'Not defined'

# Get participant's age
try:
    age = TBL[TBL['ID'] == subject]['age'].to_list()[0]
except:
    age = 'Not defined'

# Path to figures
figPath = os.path.join(os.path.dirname(out), "analysis", "vertexwise", subject, session) + "/" + subject

def report_block_template(figPath, subject='', session='', sex='', age='', featList_ctx='', featList_sctx='', featList_hipp=''):
    ctx_block = ('<p style="text-align:left;margin-left=0px;"> <a href="{figPath}_ctx-mz.png" target="_blank">'
                 '<img style="height:400px;margin-top:-100px;" src="{figPath}_ctx-mz.png"> </a> </p>' )
    sctx_block = ('<p style="text-align:left;margin-left=0px;margin-top=-500px;"> <a href="{figPath}_sctx-mz.png" target="_blank">'
                 '<img style="height:400px;margin-top:0;" src="{figPath}_sctx-mz.png">'
                 '</a></p>')
    hipp_block = (
        '<p style="text-align:left;margin-left=0px;margin-top=-500px;"> <a href="{figPath}_hipp-mz.png" target="_blank">'
        '<img style="height:400px;margin-top:0;" src="{figPath}_hipp-mz.png">'
        '</a></p>')

    # FLAIR
    if 'flair' in featList_ctx:
        flair_ctx_block = ('<p style="text-align:left;margin-left=0px;"> <a href="{figPath}_ctx-flair.png" target="_blank">'
                 '<img style="height:400px;margin-top:-100px;" src="{figPath}_ctx-flair.png"> </a> </p>' )
    else:
        flair_ctx_block = ('<p style="text-align:center">No cortical FLAIR</p>')

    if 'flair' in featList_sctx:
        flair_sctx_block = ('<p style="text-align:left;margin-left=0px;"> <a href="{figPath}_sctx-flair.png" target="_blank">'
                 '<img style="height:400px;margin-top:-100px;" src="{figPath}_sctx-flair.png"> </a> </p>' )
    else:
        flair_sctx_block = ('<p style="text-align:center">No subcortical FLAIR</p>')

    if 'flair' in featList_hipp:
        flair_hipp_block = ('<p style="text-align:left;margin-left=0px;"> <a href="{figPath}_hipp-flair.png" target="_blank">'
                 '<img style="height:400px;margin-top:-100px;" src="{figPath}_hipp-flair.png"> </a> </p>' )
    else:
        flair_hipp_block = ('<p style="text-align:center">No hippocampal FLAIR</p>')

    # qT1
    if 'qt1' in featList_ctx:
        qt1_ctx_block = ('<p style="text-align:left;margin-left=0px;"> <a href="{figPath}_ctx-qt1.png" target="_blank">'
                 '<img style="height:400px;margin-top:-100px;" src="{figPath}_ctx-qt1.png"> </a> </p>' )
    else:
        qt1_ctx_block = ('<p style="text-align:center">No cortical qT1</p>')

    if 'qt1' in featList_sctx:
        qt1_sctx_block = ('<p style="text-align:left;margin-left=0px;"> <a href="{figPath}_sctx-qt1.png" target="_blank">'
                 '<img style="height:400px;margin-top:-100px;" src="{figPath}_sctx-qt1.png"> </a> </p>' )
    else:
        qt1_sctx_block = ('<p style="text-align:center">No subcortical qT1</p>')

    if 'qt1' in featList_hipp:
        qt1_hipp_block = ('<p style="text-align:left;margin-left=0px;"> <a href="{figPath}_hipp-qt1.png" target="_blank">'
                 '<img style="height:400px;margin-top:-100px;" src="{figPath}_hipp-qt1.png"> </a> </p>' )
    else:
        qt1_hipp_block = ('<p style="text-align:center">No hippocampal qT1</p>')

    # ADC
    if 'adc' in featList_ctx:
        adc_ctx_block = ('<p style="text-align:left;margin-left=0px;"> <a href="{figPath}_ctx-adc.png" target="_blank">'
                 '<img style="height:400px;margin-top:-100px;" src="{figPath}_ctx-adc.png"> </a> </p>' )
    else:
        adc_ctx_block = ('<p style="text-align:center">No cortical ADC</p>')

    if 'adc' in featList_sctx:
        adc_sctx_block = ('<p style="text-align:left;margin-left=0px;"> <a href="{figPath}_sctx-adc.png" target="_blank">'
                 '<img style="height:400px;margin-top:-100px;" src="{figPath}_sctx-adc.png"> </a> </p>' )
    else:
        adc_sctx_block = ('<p style="text-align:center">No subcortical ADC</p>')

    if 'adc' in featList_hipp:
        adc_hipp_block = ('<p style="text-align:left;margin-left=0px;"> <a href="{figPath}_hipp-adc.png" target="_blank">'
                 '<img style="height:400px;margin-top:-100px;" src="{figPath}_hipp-adc.png"> </a> </p>' )
    else:
        adc_hipp_block = ('<p style="text-align:center">No hippocampal ADC</p>')

    # THICKNESS
    if 'thickness' in featList_ctx:
        thickness_ctx_block = ('<p style="text-align:left;margin-left=0px;"> <a href="{figPath}_ctx-thickness.png" target="_blank">'
                 '<img style="height:400px;margin-top:-100px;" src="{figPath}_ctx-thickness.png"> </a> </p>' )
    else:
        thickness_ctx_block = ('<p style="text-align:center">No cortical atrophy</p>')

    if 'thickness' in featList_sctx:
        thickness_sctx_block = ('<p style="text-align:left;margin-left=0px;"> <a href="{figPath}_sctx-thickness.png" target="_blank">'
                 '<img style="height:400px;margin-top:-100px;" src="{figPath}_sctx-thickness.png"> </a> </p>' )
    else:
        thickness_sctx_block = ('<p style="text-align:center">No subcortical atrophy</p>')

    if 'thickness' in featList_hipp:
        thickness_hipp_block = ('<p style="text-align:left;margin-left=0px;"> <a href="{figPath}_hipp-thickness.png" target="_blank">'
                 '<img style="height:400px;margin-top:-100px;" src="{figPath}_hipp-thickness.png"> </a> </p>' )
    else:
        thickness_hipp_block = ('<p style="text-align:center">No hippocampal atrophy</p>')

    report_block = (
        # =============================================================================================================#
        # =============================================== P A G E   # 1 ===============================================#
        # =============================================================================================================#
        # Title
        '<hr> <p style="margin-bottom:0;font-family:gill sans, sans-serif;text-align:center;font-size:28px;'
        'color:#181818"> Clinical report &nbsp; | &nbsp; <b> Regional changes </b> </p>'

        # Subject's ID - Session - Basic demographics
        '<p style="margin-bottom:0;margin-top:-100px;font-family:gill sans,sans-serif;text-align:center;font-size:14px;'
        'color:#181818"> <b>Subject</b>: {subject}, <b>Session</b>: {session}, <b>Sex</b>: {sex}, <b>Age</b>: {age}, '
        '<b><i>z</i>-score threshold</b>: {thr} </p>' 
        '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;'
        'color:#13365d"> <b> Blue </b> = patient < controls (note: thickness is flipped) </p>'
        '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;'
        'color:#b31b2c"> <b> Red </b> = patient > controls (note: thickness is flipped) </p> <hr>' + '<br>' +

        # Multivariate cortical changes
        '<p style="margin-bottom:-500px;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;'
        'color:#181818"> <b> Multivariate cortical changes </b> </br> (<i>{featList_ctx}</i>) </p>' + ctx_block +
        '</br>'

        '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;'
        'color:#181818"> <b> Multivariate subcortical changes </b> </br> (<i>{featList_sctx}</i>) </p>' + sctx_block + '</br>'

        '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;'
        'color:#181818"> <b> Multivariate hippocampal changes </b> </br> (<i>{featList_hipp}</i>) </p>' +
         hipp_block + '</br>'

        + '<hr>'

        # =============================================================================================================#
        # =============================================== P A G E   # 2 ===============================================#
        # ================================================= F L A I R =================================================#
        # =============================================================================================================#
        + '<p style="page-break-before:always;"></p> </br> </br> </br> </br> <hr> </br> <p style="margin-bottom:0;margin-top:0;'
          'font-family:gill sans,sans-serif;text-align:center;font-size:14px;color:#181818"> '
          '<b> Univariate cortical FLAIR </b> </br> </p>' + flair_ctx_block + '</br>'

        + '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;'
          'text-align:center;font-size:14px;color:#181818"> <b> Univariate subcortical FLAIR </b> </br> </p>'
        + flair_sctx_block + '</br>'

        + '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;'
          'text-align:center;font-size:14px;color:#181818"> <b> Univariate hippocampal FLAIR </b> </br> </p>'
        + flair_hipp_block + '</br>' + '<hr>'

        # =============================================================================================================#
        # =============================================== P A G E   # 3 ===============================================#
        # =================================================== q T 1 ===================================================#
        # =============================================================================================================#
        + '<p style="page-break-before:always;"></p> </br> </br> </br> </br> <hr> </br> <p style="margin-bottom:0;margin-top:0;'
          'font-family:gill sans,sans-serif;text-align:center;font-size:14px;color:#181818"> '
          '<b> Univariate cortical qT1 </b> </br> </p>' + qt1_ctx_block + '</br>'

        + '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;'
          'text-align:center;font-size:14px;color:#181818"> <b> Univariate subcortical qT1 </b> </br> </p>'
        + qt1_sctx_block + '</br>'

        + '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;'
          'text-align:center;font-size:14px;color:#181818"> <b> Univariate hippocampal qT1 </b> </br> </p>'
        + qt1_hipp_block + '</br>' + '<hr>'

        # =============================================================================================================#
        # =============================================== P A G E   # 4 ===============================================#
        # =================================================== A D C ===================================================#
        # =============================================================================================================#
        + '<p style="page-break-before:always;"></p> </br> </br> </br> </br> <hr> </br> <p style="margin-bottom:0;margin-top:0;'
          'font-family:gill sans,sans-serif;text-align:center;font-size:14px;color:#181818"> '
          '<b> Univariate cortical ADC </b> </br> </p>' + adc_ctx_block + '</br>'

        + '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;'
          'text-align:center;font-size:14px;color:#181818"> <b> Univariate subcortical ADC </b> </br> </p>'
        + adc_sctx_block + '</br>'

        + '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;'
          'text-align:center;font-size:14px;color:#181818"> <b> Univariate hippocampal ADC </b> </br> </p>'
        + adc_hipp_block + '</br>' + '<hr>'

        # =============================================================================================================#
        # =============================================== P A G E   # 5 ===============================================#
        # ============================================= T H I C K N E S S =============================================#
        # =============================================================================================================#
        + '<p style="page-break-before:always;"></p> </br> </br> </br> </br> <hr> </br> <p style="margin-bottom:0;margin-top:0;'
          'font-family:gill sans,sans-serif;text-align:center;font-size:14px;color:#181818"> '
          '<b> Univariate cortical atrophy </b> </br> </p>' + thickness_ctx_block + '</br>'

        + '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;'
          'text-align:center;font-size:14px;color:#181818"> <b> Univariate subcortical atrophy </b> </br> </p>'
        + thickness_sctx_block + '</br>'

        + '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;'
          'text-align:center;font-size:14px;color:#181818"> <b> Univariate hippocampal atrophy </b> </br> </p>'
        + thickness_hipp_block + '</br>' + '<hr>'

    )

    return report_block.format(figPath=figPath, subject=subject, session=session, sex=sex, age=age, thr=thr,
                               featList_ctx=featList_ctx, featList_sctx=featList_sctx, featList_hipp=featList_hipp)

# Create PDF report
static_report = ''
_static_block = report_block_template(figPath, subject=subject, session=session, sex=sex, age=age,
                                      featList_ctx=featList_ctx, featList_sctx=featList_sctx, featList_hipp=featList_hipp)
static_report += _static_block

# Utility function
def convert_html_to_pdf(source_html, output_filename):
    # open output file for writing (truncated binary)
    result_file = open(output_filename, "w+b")

    # convert HTML to PDF
    pisa_status = pisa.CreatePDF(
            source_html,                # the HTML to convert
            dest=result_file)           # file handle to recieve result

    # close output file
    result_file.close()                 # close output file

    # return True on success and False on errors
    return pisa_status.err

convert_html_to_pdf(static_report, figPath + '_regionZ_Report.pdf')

