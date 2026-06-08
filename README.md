# runeblaze.github.io

Personal site of Baqiao Liu (RuneBlaze). Built with [Astro](https://astro.build)
+ MDX, deployed to GitHub Pages.

## Stack

- **Astro** (static output) with `@astrojs/mdx`, `@astrojs/react`, `@astrojs/sitemap`
- **Iosevka Etoile** as the main Latin face (self-hosted, Latin-subset woff2 in `public/fonts/`)
- **Lost Century** palette (Kerrie Lake) — warm letterpress aesthetic, light/dark toggle

## Develop

```sh
pnpm install
pnpm dev        # http://localhost:4321
pnpm build      # → dist/
pnpm preview    # serve the build locally
```

## Structure

```
src/
  layouts/Base.astro     header nav + dark toggle, <head>, footer
  styles/global.css      the whole theme (palette, type, layout)
  pages/
    index.mdx            /          (about + résumé)
    fun.mdx              /fun/
    portfolio.mdx        /portfolio/ (works)
    404.astro
public/
  fonts/   images/   favicon.svg
```

Content is MDX, so React islands can be dropped in later for microinteractivity
(`@astrojs/react` is already wired up).

## Deploy

Pushing to `master` triggers `.github/workflows/deploy.yml`, which builds with
pnpm and publishes `dist/` via GitHub Pages.
