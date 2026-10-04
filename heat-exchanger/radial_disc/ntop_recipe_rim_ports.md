# Radial disc gyroid heat exchanger with four tangential rim ports: nTop recipe

A step-by-step process for building the model natively in nTop. All
dimensions are in mm.

## What you are building

- **Body:** a Ø124 × 36 mm disc. The disc axis is **Z** and the origin is at the disc centre.
- **Halves:** a divider plate at z = 0 splits the disc. Hot is in the upper half and cold in the lower half.
- **Ports:** four 3/8"-18 NPT male ports, all **tangential on the rim**. All four point the same way round the disc (counter-clockwise, seen from +Z), so the flow swirls in one direction.

| Port | Location | Axis direction | Height | Connects to |
|---|---|---|---|---|
| HOT IN | rim, at +X side | +Y | z = +9 | walled **duct** → hub (upper) |
| HOT OUT | rim, at +Y side | −X | z = +9 | rim ring (upper) |
| COLD IN | rim, at −Y side | +X | z = −9 | rim ring (lower) |
| COLD OUT | rim, at −X side | −Y | z = −9 | walled **duct** ← hub (lower) |

**Flow path:**
- **Hot:** HOT IN → duct → hub → **outward** through the gyroid → rim ring → HOT OUT.
- **Cold:** COLD IN → rim ring → **inward** through the gyroid → hub → duct → COLD OUT.
- Together this is radial **counter-flow**.

Block names are given as nTop's block search shows them. If your version names one differently, search for the operation (for example "cylinder", "union", "subtract", "offset").

---

## Step 0. Variables
Create each as a **Scalar (Length)** variable. To drive the model with nTop Automate later, right-click a variable and choose **Make Input**.

| Variable | Value | Meaning |
|---|---|---|
| R_OUT | 62 | disc outer radius |
| H | 36 | disc height (z = −18 … +18) |
| R_EDGE | 6 | edge rounding |
| SHELL | 2 | outer wall |
| R_HUB | 12 | hub (centre header) radius |
| R_CORE | 50 | outer radius of the gyroid core |
| R_RING | 57 | outer radius of the rim ring header |
| DIV | 2 | divider plate thickness at z = 0 |
| SEAL | 1.5 | radial width of the two seal rings |
| CELL | 6 | gyroid cell size |
| T_WALL | 0.6 | gyroid wall thickness |
| PORT_OD | 17.145 | 3/8" NPT boss diameter |
| PORT_BORE | 10 | port bore diameter |
| PORT_LEN | 14 | boss length beyond the disc |
| Z_PORT | 9 | port height (hot +9, cold −9) |
| X_T | 53.5 | tangential offset of the port axes (centre of the rim ring) |
| Y_EXIT | 31.3 | where a port axis leaves the disc = √(R_OUT² − X_T²) |
| DUCT_BORE | 8 | inner diameter of the two hub ducts |
| DUCT_WALL | 1.5 | duct wall thickness (duct OD = 11) |

---

## Step 1. Disc
1. **Cylinder**: centre (0, 0, −18), axis +Z, radius R_OUT, height H. Name it `Disc`.
2. Round the edges with R_EDGE. Use **Fillet Body**, or an **Offset Body** of −6 followed by an **Offset Body** of +6.

## Step 2. Four tangential port bosses
Each is a **Cylinder** of radius PORT_OD/2.

| Name | Start point | Axis | Length |
|---|---|---|---|
| `Boss_HotIn` | (53.5, −5.5, +9) | +Y | Y_EXIT + PORT_LEN + 5.5 |
| `Boss_HotOut` | (0, 53.5, +9) | −X | Y_EXIT + PORT_LEN |
| `Boss_ColdIn` | (0, −53.5, −9) | +X | Y_EXIT + PORT_LEN |
| `Boss_ColdOut` | (−53.5, +5.5, −9) | −Y | Y_EXIT + PORT_LEN + 5.5 |

3. **Boolean Union**: `Disc` + the four bosses, with a smoothing (blend) radius of about 3. Name it `Outer`.

## Step 3. Inner space
4. **Offset Body** `Outer` by −SHELL. Name it `Inner`.

## Step 4. Gyroid sides
5. Build the gyroid field with k = 2π/CELL = 1.0472 /mm:
   **F = sin(kx)·cos(ky) + sin(ky)·cos(kz) + sin(kz)·cos(kx)**.
   Use nTop's TPMS gyroid field with cell size CELL, or build it with math blocks from the X, Y and Z coordinate fields.
6. Use wall level **g0 ≈ 0.38**, which gives a wall of about 0.6 mm. Measure the wall in a section and adjust g0 if needed.
   - `HotSide` = **Body from Field** where F > +g0.
   - `ColdSide` = **Body from Field** where F < −g0.

## Step 5. Core, halves and seal rings
7. `CoreRing` = `Inner` ∩ annulus (radius R_HUB … R_CORE, full height).
8. `Upper` = half-space z > +DIV/2. `Lower` = half-space z < −DIV/2. Make each as a large box.
9. `Seals` = annulus (R_HUB … R_HUB+SEAL) ∪ annulus (R_CORE−SEAL … R_CORE), full height.
10. `HotCore` = (`CoreRing` ∩ `HotSide`) − (`Seals` − `Upper`).
    At the hub and rim edges, the hot channels can open only in the upper half.
11. `ColdCore` = (`CoreRing` ∩ `ColdSide`) − (`Seals` − `Lower`).

## Step 6. Headers
12. `Hub` = `Inner` ∩ **Cylinder** (radius R_HUB+0.5, full height).
13. `Ring` = `Inner` ∩ annulus (R_CORE−0.5 … R_RING, full height).
14. `HotVoid` = `HotCore` ∪ (`Hub` ∩ `Upper`) ∪ (`Ring` ∩ `Upper`).
15. `ColdVoid` = `ColdCore` ∪ (`Hub` ∩ `Lower`) ∪ (`Ring` ∩ `Lower`).

## Step 7. Hub ducts (the key step for rim-only ports)
HOT IN and COLD OUT have to reach the hub. Each uses an L-shaped, walled duct.

**Duct walls** (OD = DUCT_BORE + 2·DUCT_WALL = 11):

16. `DuctWall_Hot` = **Boolean Union** of:
    - a radial leg: **Cylinder** from (R_HUB−1, 0, +9) to (59, 0, +9), axis +X, radius 5.5;
    - a tangential leg: **Cylinder** from (53.5, −5.5, +9), axis +Y, length Y_EXIT + 5.5, radius 5.5.
17. `DuctWall_Cold` is the same shape rotated 180° about Z and moved to z = −9:
    - radial leg along −X;
    - tangential leg from (−53.5, +5.5, −9) along −Y.

**Cut the walls out of both fluids.** This step is what keeps the duct sealed from the core and from the rim ring it crosses:

18. `HotVoid` = `HotVoid` − (`DuctWall_Hot` ∪ `DuctWall_Cold`).
19. `ColdVoid` = `ColdVoid` − (`DuctWall_Hot` ∪ `DuctWall_Cold`).

**Duct bores:**

20. `DuctBore_Hot` = **Boolean Union** of:
    - a radial bore, radius 4, from (R_HUB−3, 0, +9) to (53.5, 0, +9). It opens into the hub;
    - a tangential bore, radius 4, from (53.5, 0, +9) along +Y to the boss tip;
    - the port bore, radius PORT_BORE/2, from y = Y_EXIT−1 to the boss tip.
21. `DuctBore_Cold` is the same, mirrored, at z = −9.
22. `HotVoid` = `HotVoid` ∪ `DuctBore_Hot`. `ColdVoid` = `ColdVoid` ∪ `DuctBore_Cold`.

## Step 8. Ring ports (HOT OUT, COLD IN)
23. `Bore_HotOut` = **Cylinder** radius PORT_BORE/2, inside `Boss_HotOut` from its start point to its tip, **minus** a **Cylinder** of radius R_RING−1 (full height).
24. `Bore_ColdIn` = the same, inside `Boss_ColdIn`.

    This clip matters. The Ø10 bore is wider than the 7 mm ring, so without it the bore cuts through the seal ring and hot touches cold.
25. `HotVoid` = `HotVoid` ∪ `Bore_HotOut`. `ColdVoid` = `ColdVoid` ∪ `Bore_ColdIn`.

## Step 9. Final part
26. `HotVoid` = `HotVoid` ∩ `Outer`. `ColdVoid` = `ColdVoid` ∩ `Outer`.
27. `HX` = **Boolean Subtract**: `Outer` − (`HotVoid` ∪ `ColdVoid`).

## Step 10. Checks in nTop (do these before printing)
- **Leak test:** **Boolean Intersect** `HotVoid` ∩ `ColdVoid`. The result must be **empty**. Then repeat it with each body offset by +0.1 mm. That second test catches walls of zero thickness, which a plain intersect misses.
- **Section at z = +9:** red path HOT IN → L-duct → hub, gyroid, rim ring → HOT OUT. The duct must have a grey wall all round, including where it crosses the rim ring.
- **Section at z = −9:** the mirror image for cold.
- **Section at y = 0:** check the divider plate at z = 0, the hub split into upper and lower parts, and the seal rings at the hub and rim.
- **Wall thickness:** measure the gyroid wall in a section. It should be about T_WALL; adjust g0 if not.
- **Fluid volumes:** `HotVoid` and `ColdVoid` should be about equal, roughly 130 cm³ each.

## Step 11. Threads
I don't know of a dedicated NPT thread block in nTop. Two options:
- Leave the bosses as plain Ø17.145 cylinders and cut the 3/8"-18 NPT thread with a die after printing, or
- Union in a thread body modelled in CAD.

## Design notes
- **Counter-flow:** hot spreads outward while cold converges inward through the same gyroid.
- **Velocity toward the rim:** the flow area grows about 4× from hub to rim, so velocity, and with it the heat-transfer coefficient, falls toward the rim. A finer cell toward the rim (graded CELL) would even this out.
- **Ducts:** each duct takes a narrow strip of core at its height. The gyroid flows around it.
- **Swirl:** all four tangential ports turn the same way. Swirl in the rim ring helps spread the flow round the full circumference.
- **Pre-check:** before writing this recipe, I checked this exact layout in a private test model with a 0.25 mm grid:
  - hot and cold touch nowhere;
  - each fluid is one connected volume of about 130 cm³;
  - all four ports open into the correct fluid.

  I have **not** built it in nTop itself, so do the Step 10 checks there.
