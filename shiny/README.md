# Stratified Sampling Simulation Explorer

A small Shiny app to explore `med_all_res.RData` (simulation results comparing
stratification methods across datasets, strata counts, and seeds).

## What's in this folder

- `app.R` — the Shiny app (single-file: UI + server)
- `med_all_res.RData` — the data (must contain a data frame named `res`)

## Run it locally

```r
install.packages(c("shiny", "dplyr", "ggplot2", "DT", "tidyr"))
shiny::runApp("app.R")
```

## Put it on GitHub (just the code)

```bash
cd shiny_app
git init
git add .
git commit -m "Initial Shiny app for simulation results"
git branch -M main
git remote add origin https://github.com/<your-username>/<repo-name>.git
git push -u origin main
```

This makes your code and data visible/downloadable on GitHub, but **GitHub
itself does not execute R code**, so the app won't be "live" from a repo
alone. To let other people actually click around in the app in a browser,
pick one of the two options below.

---

## Option A — shinyapps.io (real R backend, easiest, free tier)

This runs your actual R code on Posit's servers.

1. Create a free account at https://www.shinyapps.io/
2. In R:
   ```r
   install.packages("rsconnect")
   rsconnect::setAccountInfo(name   = "<your account name>",
                             token  = "<token from shinyapps.io dashboard>",
                             secret = "<secret from shinyapps.io dashboard>")
   rsconnect::deployApp("path/to/shiny_app")
   ```
3. You'll get a public URL like `https://<you>.shinyapps.io/shiny_app/` to share.
4. You can keep the GitHub repo as the source of truth and redeploy after each
   `git push` by re-running `rsconnect::deployApp()`.

I'm confident this general workflow is correct, but exact free-tier limits
(active hours/month, number of apps) change over time — check current
numbers at https://www.shinyapps.io/ before planning around them.

## Option B — shinylive on GitHub Pages (no server, runs in the browser)

`shinylive` compiles your Shiny app to WebAssembly so it runs client-side —
this is the way to get something genuinely hosted "on git" (GitHub Pages)
with no separate server. I'm less certain of the exact current CLI/package
names since this tooling has evolved — please cross-check against
https://posit-dev.github.io/r-shinylive/ before running these commands.

```r
install.packages("shinylive")
shinylive::export("shiny_app", "docs")  # "docs" folder for GitHub Pages
```

Then push the `docs/` folder to your repo and enable GitHub Pages
(Settings → Pages → source: `main` branch, `/docs` folder).

Caveats worth knowing about shinylive before you commit to it:
- Every visitor downloads a WebAssembly R runtime in their browser (larger
  initial load, slower first render).
- Not all R packages work in this environment — `dplyr`, `ggplot2`, `DT`,
  and `tidyr` are generally supported, but if you add anything unusual,
  verify it's shinylive-compatible first.
- Large `.RData` files will be shipped to every visitor's browser as a
  static asset — fine for small result tables like this one (a few KB),
  but not a good fit if your real dataset grows into the tens of MB.

## Which to choose

- Want the simplest path and don't mind a small external dependency →
  **shinyapps.io**.
- Want everything, including hosting, to live in one GitHub repo with no
  external account →**shinylive + GitHub Pages**, given the caveats above.
