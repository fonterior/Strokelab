# Third-Party Notices

StrokeLab uses third-party Python packages. Each package remains subject to its own license. The StrokeLab license does not replace, modify, or relicense these dependencies.

This file lists the direct dependencies declared in `requirements.txt` and `requirements-optional.txt`. It is not a complete software bill of materials (SBOM): package managers may install additional transitive dependencies, and the exact dependency tree can change when package versions change.

The StrokeLab repository does not intentionally vendor the packages listed below. They are installed separately by the user with `pip`. If a binary release bundles any of these packages, include the applicable license texts and notices with that release and review the package's current distribution terms.

## Required dependencies

| Package | Purpose in StrokeLab | Upstream project | License information |
|---|---|---|---|
| `fonttools` | Reading and processing font files and vector outlines | <https://github.com/fonttools/fonttools> | MIT; see the upstream [`LICENSE`](https://github.com/fonttools/fonttools/blob/main/LICENSE). |
| `PySide6` | Qt for Python bindings used by the graphical interface | <https://wiki.qt.io/Qt_for_Python> | The current PyPI package metadata declares LGPL-3.0-only, GPL-2.0-only, or GPL-3.0-only options; Qt also offers a separate commercial license. See the [Qt for Python licensing overview](https://doc.qt.io/qtforpython-6/) and [Qt's LGPL/GPL obligations](https://www.qt.io/development/open-source-lgpl-obligations). |

PySide6 is a special case. The PySide6 6.12.0 package metadata lists LGPL-3.0-only, GPL-2.0-only, or GPL-3.0-only licensing options; Qt also offers a separate commercial license. The exact terms should be checked against the version and modules actually installed. Some Qt modules may have different availability or licensing. Before distributing an executable that bundles PySide6/Qt, identify the exact Qt modules included and verify that the selected license permits the intended distribution. Under LGPL, obligations can include providing the complete corresponding source code for the library, prominent license notices, and a way for recipients to replace/relink the library and run the resulting application. Qt's guidance also states that recipients' LGPL freedoms must not be restricted. This notice is an orientation, not a legal assessment of a particular binary build.

## Optional dependencies

| Package | Purpose in StrokeLab | Upstream project | License information |
|---|---|---|---|
| `skia-pathops` | Optional boolean/path operations, including overlap removal | <https://github.com/fonttools/skia-pathops> | BSD 3-Clause; see the upstream [`LICENSE`](https://github.com/fonttools/skia-pathops/blob/main/LICENSE). |
| `glyphsLib` | Optional Glyphs source/export support | <https://github.com/googlefonts/glyphsLib> | Apache License 2.0; see the upstream [`LICENSE`](https://github.com/googlefonts/glyphsLib/blob/main/LICENSE). |
| `ufoLib2` | Optional UFO data model and export support | <https://github.com/fonttools/ufoLib2> | Apache License 2.0; see the upstream [`LICENSE`](https://github.com/fonttools/ufoLib2/blob/main/LICENSE). |

## Important distinction between licenses

The PolyForm Noncommercial License applies to the StrokeLab software authored for this repository. It does not change the licenses of FontTools, PySide6/Qt, skia-pathops, glyphsLib, ufoLib2, or their transitive dependencies.

The permissive licenses used by some dependencies may allow commercial use of those dependency packages. That permission does not grant commercial-use rights to StrokeLab itself. Commercial use of StrokeLab requires separate permission from its author.

## Changes and release checks

When changing dependencies or preparing a release:

1. Check the license metadata and license files for the exact package versions being shipped.
2. Inspect transitive dependencies, not only the five direct dependencies listed here.
3. Include required third-party license texts and notices in any release that bundles third-party components.
4. For PySide6/Qt binary bundles, perform a separate LGPL/GPL compliance review before distribution.
