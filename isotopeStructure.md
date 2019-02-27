## Notes and variables for isotope data
Each isotope contains the following named variables, and can be called using the $ syntax in R.  For example, `Isotopes$AC225$t12` returns `10` as a numeric value

- **isotope**: Character string as needed for searching (all caps)
- **decayLevel**: The number of decays from the initial isotope
- **A:** Mass number
- **Z:** Atomic number
- **symb:** Element symbol
- **masterYield:** The (fractional) quantity of the isotope relative to the initial isotope
- **t12:** The isotope's half-life in units of days
- **SA** The isotope's specific activity in units of mCi / ug

**Decays:** A list with length corresponding to the number of decay branches originating from the isotope. The names of the list items within **Decays** identifies the nature of the decay as sourced from the JEFF data file. Current possibilities are:
- **Alpha**
- **Beta**
- **Positron**
- **EC** (electron capture)
- **IT** (internal transition)

Within each sublist (e.g. **Alpha**), the following properties can be found:
- **Q:** The decay energy in keV
- **branchYield:** The fraction of decays that follow the containing branches
- **BrQ:** The decay energy Q adjusted for branch yield, also in keV. This parameter is calculated by `Q * branchYield`
- **daughter:** The identity of the daughter product resulting from the containing decay mechanism. Note that the lower-case symbol is not compatible with direct search of the data file (use `toupper`)

For example, `Isotopes$AC225$Decays$Alpha$daughter` will return `"221Fr"` as a text string.
