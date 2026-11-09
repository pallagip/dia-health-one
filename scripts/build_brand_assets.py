from pathlib import Path
import subprocess
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"assets/brand"
OUT.mkdir(parents=True,exist_ok=True)
V="#7155bd"; V2="#927dd4"; L="#c8e5a8"; I="#201f27"; IV="#f8f8f5"; MUT="#686671"
def put(name,svg,export=None):
 p=OUT/name
 p.write_text(svg,encoding="utf-8")
 for png,w,h in export or []:
  subprocess.run(["rsvg-convert","-w",str(w),"-h",str(h),"-o",str(OUT/png),str(p)],check=True)
  assert Image.open(OUT/png).size==(w,h)
  print(png,w,h,(OUT/png).stat().st_size)
def mark(x,y,size):
 scale=size/512
 return f'<g transform="translate({x} {y}) scale({scale})"><rect width="512" height="512" rx="132" fill="{V}"/><path d="M138 112h105c91 0 145 56 145 144S334 400 243 400H138V112Zm76 67v154h29c45 0 70-27 70-77s-25-77-70-77h-29Z" fill="#fff"/><circle cx="374" cy="138" r="28" fill="{L}"/></g>'
def svg(w,h,inner,defs=""):
 return f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" role="img"><defs>{defs}</defs>{inner}</svg>'
font='font-family="Manrope,Liberation Sans,Arial,sans-serif"'
mono='font-family="DM Mono,DejaVu Sans Mono,monospace"'
base=(ROOT/"assets/dia-health-one-mark.svg").read_text()
put("dia-mark-master.svg",base,[("linkedin-logo-400.png",400,400),("dia-logo-1024.png",1024,1024)])
word=mark(0,7,156)+f'<text x="185" y="118" fill="{I}" font-size="86" font-weight="800" letter-spacing="-3" {font}>DIA</text><text x="366" y="118" fill="{V}" font-size="86" font-weight="800" letter-spacing="-3" {font}>HEALTH</text><text x="757" y="118" fill="{I}" font-size="86" font-weight="800" letter-spacing="-3" {font}>ONE</text>'
put("wordmark-dark.svg",svg(1200,180,word),[("wordmark-dark.png",1200,180)])
word2=mark(0,7,156)+f'<text x="185" y="118" fill="white" font-size="86" font-weight="800" letter-spacing="-3" {font}>DIA</text><text x="366" y="118" fill="{L}" font-size="86" font-weight="800" letter-spacing="-3" {font}>HEALTH</text><text x="757" y="118" fill="white" font-size="86" font-weight="800" letter-spacing="-3" {font}>ONE</text>'
put("wordmark-light.svg",svg(1200,180,word2))
defs=f'<linearGradient id="bg" x1="0" y1="0" x2="1" y2="1"><stop stop-color="#201f27"/><stop offset=".64" stop-color="#3b2d62"/><stop offset="1" stop-color="#6651a5"/></linearGradient><radialGradient id="glow"><stop stop-color="{V2}" stop-opacity=".55"/><stop offset="1" stop-color="{V2}" stop-opacity="0"/></radialGradient><pattern id="dots" width="68" height="68" patternUnits="userSpaceOnUse"><circle cx="3" cy="3" r="2.4" fill="#cec2f9" opacity=".35"/></pattern>'
parts=[f'<rect width="4200" height="700" fill="url(#bg)"/>','<circle cx="3570" cy="310" r="650" fill="url(#glow)"/>','<rect x="2880" width="1320" height="700" fill="url(#dots)" opacity=".8"/>']
for r in (180,310,450,600,780):
 parts.append(f'<circle cx="3630" cy="370" r="{r}" fill="none" stroke="#e3dff5" stroke-width="3" opacity=".18"/>')
parts.extend([
 f'<path d="M2530 615 C2910 420,3030 520,3340 350 S3800 230,4200 70" fill="none" stroke="{L}" stroke-width="7" opacity=".7"/>',
 f'<path d="M2450 690 C2940 480,3270 750,3670 490 S4030 395,4300 310" fill="none" stroke="{V2}" stroke-width="24" opacity=".45"/>',
 f'<circle cx="3340" cy="350" r="21" fill="{L}"/>',f'<circle cx="3825" cy="224" r="21" fill="{V2}"/>',mark(410,254,180),
 f'<text x="660" y="335" fill="white" font-size="165" font-weight="800" letter-spacing="-7" {font}>DIA HEALTH ONE</text>',
 f'<text x="669" y="425" fill="{L}" font-size="53" font-weight="600" letter-spacing="1" {font}>Computational insulin discovery</text>',
 f'<text x="670" y="498" fill="#d3cde3" font-size="31" letter-spacing="4" {mono}>EVIDENCE-LED  /  NON-CLINICAL RESEARCH</text>',
 '<path d="M670 546 H1920" stroke="#ffffff" opacity=".2" stroke-width="3"/>'
])
put("linkedin-cover-master.svg",svg(4200,700,"".join(parts),defs),[("linkedin-cover-4200x700.png",4200,700),("linkedin-cover-1512x256.png",1512,256),("linkedin-cover-1128x191.png",1128,191)])
parts=[f'<rect width="1200" height="627" fill="{IV}"/>',f'<circle cx="1130" cy="180" r="400" fill="#e8e1f5"/>',f'<circle cx="1090" cy="175" r="245" fill="none" stroke="#d5caea" stroke-width="2"/>',f'<circle cx="1090" cy="175" r="169" fill="none" stroke="#d5caea" stroke-width="2"/>',f'<circle cx="1010" cy="483" r="140" fill="{L}" opacity=".55"/>',f'<path d="M835 528 C925 423,971 493,1033 405 S1140 388,1220 301" fill="none" stroke="{V}" stroke-width="7" opacity=".7"/>',f'<circle cx="1033" cy="405" r="16" fill="{V}"/>',mark(72,60,74),f'<text x="167" y="112" fill="{I}" font-size="30" font-weight="800" letter-spacing="1" {font}>DIA HEALTH ONE</text>',f'<text x="72" y="254" fill="{I}" font-size="70" font-weight="800" letter-spacing="-3" {font}>RESEARCH, WITH</text>',f'<text x="72" y="341" fill="{V}" font-size="70" font-weight="800" letter-spacing="-3" {font}>EVIDENCE &amp; LIMITS.</text>',f'<text x="77" y="416" fill="{MUT}" font-size="26" {font}>Exploring glucose-responsive insulin concepts</text>',f'<rect x="72" y="479" width="235" height="54" rx="27" fill="{L}"/>',f'<text x="190" y="513" text-anchor="middle" fill="{I}" font-size="17" font-weight="800" letter-spacing="1.5" {mono}>NON-CLINICAL</text>',f'<line x1="72" y1="564" x2="1128" y2="564" stroke="#dedbe1" stroke-width="2"/>',f'<text x="75" y="596" fill="{MUT}" font-size="18" {mono}>DIAHEALTH.ONE</text>',f'<text x="1122" y="596" text-anchor="end" fill="{MUT}" font-size="17" {mono}>RESEARCH UPDATE / 01</text>']
put("linkedin-news-1200x627.svg",svg(1200,627,"".join(parts)),[("linkedin-news-1200x627.png",1200,627)])
# Keep the existing website OpenGraph graphic in the same design system.
site_og=ROOT/"assets/dia-health-one-social.svg"
site_og.write_text(svg(1200,630,"".join(parts)),encoding="utf-8")
subprocess.run(["rsvg-convert","-w","1200","-h","630","-o",str(ROOT/"assets/dia-health-one-social.png"),str(site_og)],check=True)
sq=[f'<rect width="1080" height="1080" fill="{IV}"/>',f'<rect width="1080" height="32" fill="{V}"/>',f'<circle cx="1010" cy="180" r="270" fill="#e7e0f6"/>',f'<circle cx="1010" cy="180" r="200" stroke="#cdbfeb" stroke-width="2" fill="none"/>',mark(64,83,92),f'<text x="182" y="147" fill="{I}" font-size="41" font-weight="800" {font}>DIA HEALTH ONE</text>',f'<text x="71" y="291" fill="{V}" font-size="26" font-weight="800" letter-spacing="5" {mono}>RESEARCH NOTE</text>',f'<text x="65" y="415" fill="{I}" font-size="90" font-weight="800" letter-spacing="-4" {font}>Research update.</text>',f'<text x="70" y="485" fill="{MUT}" font-size="30" {font}>Evidence first. Limitations explicit.</text>',f'<rect x="64" y="562" width="950" height="324" rx="30" fill="#fff" stroke="#dfd8eb" stroke-width="3"/>',f'<circle cx="117" cy="631" r="11" fill="{V}"/><text x="150" y="641" fill="{I}" font-size="32" font-weight="700" {font}>What we tested</text>',f'<line x1="115" y1="690" x2="955" y2="690" stroke="#ece9ee" stroke-width="2"/>',f'<circle cx="117" cy="758" r="11" fill="{V2}"/><text x="150" y="768" fill="{I}" font-size="32" font-weight="700" {font}>What we observed</text>',f'<circle cx="117" cy="842" r="11" fill="#86b267"/><text x="150" y="852" fill="{I}" font-size="32" font-weight="700" {font}>What remains unknown</text>',f'<rect x="64" y="928" width="325" height="62" rx="31" fill="{L}"/>',f'<text x="226" y="968" text-anchor="middle" fill="{I}" font-size="22" font-weight="800" letter-spacing="2" {mono}>NON-CLINICAL</text>',f'<text x="1010" y="973" fill="{MUT}" text-anchor="end" font-size="22" {mono}>DIAHEALTH.ONE</text>']
put("linkedin-research-template-1080.svg",svg(1080,1080,"".join(sq)),[("linkedin-research-template-1080.png",1080,1080)])
