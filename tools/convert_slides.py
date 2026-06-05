"""Convert multiplierless-FIR-slides.md to Remark.js HTML."""
import re

with open("multiplierless-FIR-slides.md", "r", encoding="utf-8") as f:
    text = f.read()

# 1. Strip YAML front matter
text = re.sub(r"^---\n.*?\n---\n", "", text, flags=re.DOTALL)

# 2. Strip HTML comments
text = re.sub(r"<!--.*?-->\n?", "", text)

# 3. Split into slides
slides = [s.strip() for s in text.split("\n---\n") if s.strip()]

# Extract layout directive (must be at top before all slides)
layout = ""
if slides and slides[0].startswith("layout: true"):
    layout = slides[0] + "\n\n---\n\n"
    slides = slides[1:]

# 4. Annotate slides
annotated = []
for i, s in enumerate(slides):
    head = "\n".join(s.split("\n")[:3])  # check first 3 lines for existing class
    has_class = any(c in head for c in ("nord-dark", "nord-light"))
    if i == 0 and not has_class:
        annotated.append(f"count: false\nclass: nord-dark, middle, center\n\n{s}")
    else:
        annotated.append(s)

result = layout + "\n\n---\n\n".join(annotated)

# 5. Wrap in HTML
HTML = r"""<!doctype html>
<html>
  <head>
    <title>Multiplierless FIR Filter Design Toolchain</title>
    <meta charset="utf-8" />
    <meta name="viewport"
      content="user-scalable=no,initial-scale=1,maximum-scale=1,minimum-scale=1,width=device-width" />
    <link rel="stylesheet" type="text/css" href="../katex/katex.min.css" />
    <link rel="stylesheet" type="text/css" href="../css/spaces.css" />
    <link rel="stylesheet" type="text/css" href="../css/slides.css" />
    <link rel="stylesheet" type="text/css" href="../css/nord-dark.css" />
    <link rel="stylesheet" type="text/css" href="../css/nord-light.css" />
    <link rel="stylesheet" type="text/css" href="../css/font-nord.css" />
    <link rel="stylesheet" type="text/css" href="../css/bg-nord.css" />
    <link rel="stylesheet" type="text/css" href="../css/style.css" />
  </head>
  <body>
    <textarea id="source">
{slides}
    </textarea>
    <script src="../js/remark.min.js"></script>
    <script src="../katex/katex.min.js" type="text/javascript"></script>
    <script src="../katex/contrib/auto-render.min.js" type="text/javascript"></script>
    <script src="../js/mermaid.min.js"></script>
    <script type="text/javascript">
      renderMathInElement(document.getElementById('source'), {
        delimiters: [
          { left: '$$', right: '$$', display: true },
          { left: '$', right: '$', display: false },
        ],
      });
      var slideshow = remark.create({
        ratio: '16:10',
        highlightStyle: 'tomorrow-night-blue',
        highlightLines: true,
        countIncrementalSlides: false,
        navigation: {
          scroll: false,
          touch: true,
          click: false,
        },
      });
    </script>
    <script src="../js/mermaid-init.js"></script>
  </body>
</html>"""

output = HTML.replace("{slides}", result)

outpath = "../../luk036.github.io/cvx/multiplierless-FIR-remark.html"
with open(outpath, "w", encoding="utf-8") as f:
    f.write(output)

print(f"Done! {len(slides)} slides → {outpath}")
