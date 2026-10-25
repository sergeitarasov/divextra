"""Split the compiled vector artwork and crop to ink bounds plus 3 pt padding."""
from pathlib import Path
import re
import subprocess
from pypdf import PdfReader, PdfWriter
from pypdf.generic import RectangleObject

folder = Path(__file__).resolve().parent
source = folder / "icon-pages.pdf"
names = ["change-0-to-1", "change-1-to-0", "split-upper-0", "split-upper-1",
         "split-lower-0", "split-lower-1", "extinction", "nothing"]
result = subprocess.run(["gs", "-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=bbox",
                         str(source)], capture_output=True, text=True, check=True)
boxes = re.findall(r"%%HiResBoundingBox: (.*)", result.stderr)
reader = PdfReader(source)
assert len(reader.pages) == len(boxes) == len(names)
(folder / "icons").mkdir(exist_ok=True)
for page, box, name in zip(reader.pages, boxes, names):
    left, bottom, right, top = map(float, box.split())
    bounds = RectangleObject([left - 3, bottom - 3, right + 3, top + 3])
    page.mediabox = bounds
    page.cropbox = bounds
    writer = PdfWriter()
    writer.add_page(page)
    with (folder / "icons" / f"{name}.pdf").open("wb") as output:
        writer.write(output)
print("Exported 8 tightly cropped vector PDFs.")
