# Cylindrical gyroid heat exchanger: VS Code → nTop

## Files

| File | Purpose |
|---|---|
| `cyl_gyroid_hx.py` | Generates the model and writes `HX_solid.stl`, `Hot_fluid.stl` and `Cold_fluid.stl` |
| `run_ntop_automate.py` | Optional: runs the generator, then runs an nTop notebook headless (nTop Automate / nTopCL) on each mesh |
| `cyl_hx_concept.py` | Concept renders only (coarse, thick walls) |
| `requirements.txt` | Python packages |

## 1. Set up VS Code

1. Install Python 3.10 or newer and the VS Code *Python* extension.
2. Open this folder in VS Code and open a terminal:
   ```
   python -m venv .venv
   .venv\Scripts\activate          # Windows  (macOS/Linux: source .venv/bin/activate)
   pip install -r requirements.txt
   ```

## 2. Generate the model

```
python cyl_gyroid_hx.py --preview              # colour sections -> hx_preview.png (seconds)
python cyl_gyroid_hx.py --voxel 0.4 --check    # quick coarse test
python cyl_gyroid_hx.py                        # full export, voxel 0.2 mm
```

All dimensions are in the PARAMETERS block at the top of `cyl_gyroid_hx.py`. The ones you will most likely change:

| Parameter | Default | Meaning |
|---|---|---|
| `CELL` | 8 × 8 × 8 mm | Gyroid cell size |
| `T_WALL` | 0.6 mm | Gyroid wall thickness. Must be at least about 3 × `VOXEL`, otherwise the walls break up in the mesh |
| `VOXEL` | 0.2 mm | Mesh resolution. 0.15 mm gives crisper threads |
| `R`, `HB` | 70, 90 mm | Cylinder radius and body height |
| `Z_APEX` | 20 mm | Arch apex height inside the body |
| `LEG_H`, `SKIRT` | 40, 28 mm | Tip depth and shoulder depth; the thread length is the difference (12 mm) |
| `BORE_R` | 5 mm | Port bore radius |

A full run at 0.2 mm is about 400 million grid points. It is processed in x-slices (`CHUNK`), so memory stays moderate, but expect tens of minutes and STL files of several GB. If you run out of RAM, lower `CHUNK`.

## 3. Bring it into nTop (manual)

1. **Import Mesh**: select `hx_output/HX_solid.stl`.
2. **Implicit Body from Mesh**: use that mesh as input.
3. Use the implicit body like any nTop body: section it, analyse it, or mesh it.
4. For printing, use **Mesh from Implicit Body**, then export (STL or 3MF), or use the imported mesh directly.
5. For CFD, import `Hot_fluid.stl` and `Cold_fluid.stl` the same way. They are the two fluid domains.

## 4. Notebook for nTop Automate (optional, one-time GUI setup)

nTop Automate runs an existing notebook whose variables are marked as inputs. It cannot create a notebook from scratch.

1. Create a new notebook.
2. Add a *Text* or *File Path* variable named **`Mesh Path`**.
3. Add a *Text* variable named **`Output Path`**.
4. Make both variables inputs (right-click → *Make Input*).
5. Add **Import Mesh** (from `Mesh Path`) → **Implicit Body from Mesh** → **Mesh from Implicit Body** → an export block writing to `Output Path`.
6. Save it as `hx_import.ntop`.
7. Run from VS Code:
   ```
   python run_ntop_automate.py --ntop "C:/Program Files/nTopology/nTopology/nTopCL.exe" --notebook hx_import.ntop
   ```

The JSON input format can vary between nTop versions. If nTopCL rejects the default format, generate the input template for your notebook with nTopCL (`nTopCL.exe --help` lists the option) and pass it with `--template your_template.json`.

Reference: [Running nTop Automate in Python scripts](https://support.ntop.com/hc/en-us/articles/360052703693-Running-nTop-Automate-in-Python-scripts), [Preparing an nTop Notebook for nTop Automate](https://support.ntop.com/hc/en-us/articles/360052833053-Preparing-an-nTop-Notebook-for-nTop-Automate).

## Design summary

- **Body:** a Ø140 mm cylinder from z = −28 to z = 90 mm, top edge rounded R10, 2 mm shell.
- **Arch and pillars:** a groin-vault arch (two crossing elliptical arches, apex at z = 20 mm) leaves four pillars at 45°, 135°, 225° and 315°. They taper slowly to Ø22 mm at z = −28 mm.
- **Ports:** 3/8"-18 NPT male thread (ASME B1.20.1 basic profile, 12 mm long) to the tip at z = −40 mm, with a Ø10 mm bore.
  - **Hot:** 45° and 225° pillars.
  - **Cold:** 135° and 315° pillars.
- **Headers:** inside the pillars, following the arch with a 2 mm wall, split into quadrants.
- **Core:** a sheet gyroid from z = 30 mm to the top. A seal slab at its base lets each fluid's labyrinth open only into its own two quadrants.

## Known limits

- **Gyroid wall thickness:** the wall is built from the nodal approximation divided by its gradient, so its thickness is close to `T_WALL` but not exact.
- **Threads:** the thread is a basic profile from implicit geometry. Check the printed fit with an NPT ring gauge or a fitting.
- **Thin thread wall:** the wall between the thread root and the bore is about 2.2 mm. Use `BORE_R = 4.5` for about 2.7 mm.
- **Flow not simulated:** flow distribution and pressure drop have not been checked with CFD.
- **Trapped micro-pockets:** a coarse test shows a few single-voxel fluid pockets (about 0.01 mm³) trapped at the shell. They are harmless.
