
## 2D Graphics Libraries

### TText

- The character position was not correct with the Cocoa backend.
  (see https://sft.its.cern.ch/jira/browse/ROOT-6561)
- Interactive update of `TText` position did not work in NDC mode.
  (se https://sft.its.cern.ch/jira/browse/ROOT-7284).

### TLegend

- Use the new `TStyle` global attribute `gStyle->GetLegendTextSize()` to set the
  legend item text size. If this value is 0 and if the text size directly set on
  the `TLegend` object is also 0, then the text size is automatically computed to
  fit the legend box. If `gStyle->GetLegendTextSize()` is non equal to 0 and if
  text size  directly set on the `TLegend` object is 0, then the `gStyle` value is
  used to draw the legend text. If the text size directly set on the `TLegend`
  object is not null, then it is used to draw the legend text.

### TTexDump

- The hollow fill style was not rendered correctly.
  (see https://sft.its.cern.ch/jira/browse/ROOT-6841)
- Better line width matching with screen and pdf output.
- Text color was ignored. It was always black.
- Text color was ignored. It was always black.
- The underscore `_` produced an error outside the TeX math context.
- Fix an issue with transparent pads.
- Implement transparent colors using TiKZ "opacity".
- Implement `TStyle::SetLineScalePS()` to control le global basic line width.
- Offer 0 as line width option. Useful to make a line invisible.

### TPostScript

- Small fix for fill patterns 1, 2 and 3.
- With `TMathtext`, only the fonts really used are now loaded in the PostScript
  file. Typically it reduces the file size by a factor 10 (compare to the previous
  implementation) for normal plots with math formulae and greek characters.
- Offer 0 as line width option. Useful to make a line invisible.

### TPDF

- When a text size was equal or smaller than 0 the PDF file was corrupted.
- Small fix for fill patterns 1, 2 and 3.
- When printing a coloured 2D histograms (with option COLZ) into a PDF or PostScript
  file, the preview on screen using many standard PDF previewer tools showed very
  thin white lines between the bins as well as in the color palette.
  This made very ugly the final output.
  This problem is due to bad implementation of anti-aliasing in these previewers.
  A way to bypass this issue was to turn off the anti-aliasing in the previewer
  but then the rest of the document does not look nice. This problem is now bypassed
  with a fix in both PDF and PostScript output.
- Offer 0 as line width option. Useful to make a line invisible.

### TSVG

- Use float numbers instead of integer to describe graphics paths to avoid
  rounding problems.
- Implement missing math symbols.
- Offer 0 as line width option. Useful to make a line invisible.

### TASImage

- In the animated gif it is now possible to specify the delay between the last
  image and the fist image in case of infinite loop. (Fernando Hueso González
  <f.gonzalez@hzdr.de>).

### TPadPainter
- Offer 0 as line width option. Useful to make a line invisible.

### TPad

- In `TPad::ShowGuidelines` the number of guide lines is limited to 15. Above
  that they become useless.
- Print a warning if one of the pad limit is a NaN.

### TCanvas

- Make sure that "/" and "." are not part of the method name when a canvas is
 saved as a .C file.

### TLatex

- With the Cocoa backend the PDF and PS output produced miss-aligned exponents
  because the `GetTextExtend` method behaved differently in batch mode and "screen"
  mode. This is now fixed. See http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=18883
- Improve the square-root drawing in case it is small.
- Better adjustment of the tilde accent position in case of Cocoa backend.

### TMathText

- `\mu` is now working for Postscript output.
- `\splitline` is now implemented.

### Cocoa backend

- Line width and line style were not applied on boxes.

