# cellgen website

Static project page for the cellgen simulated 4D cell microscopy dataset. No build step.

```
index.html      the page (content lives in the FAMILIES array at the bottom)
videos/         33 MIP movies, 3 per system (raw | clean)
posters/        one still per movie
tools/          mip_video_bleachcomp.py, used to render the 9 viral-CPE movies
```

## Preview locally

```bash
python -m http.server 8000      # then open http://localhost:8000
```

## Publish with GitHub Pages

```bash
git init -b main && git add -A && git commit -m "cellgen project page"
git remote add origin git@github.com:<user>/<repo>.git
git push -u origin main
```

Then on GitHub: Settings → Pages → Source "Deploy from a branch" → `main` / `/ (root)`.
The page appears at `https://<user>.github.io/<repo>/`.

## Add or swap a video

1. Copy the MP4 (H.264, yuv420p) to `videos/<name>.mp4` and a still to `posters/<name>.jpg`.
2. Edit the matching `V("<name>", "<label>", "<caption>")` entry in `index.html`.

The eight shorter systems use the movies from `gen/mip_video.py` unchanged (three SARS-CoV-2 movies
were re-encoded at CRF 27 to stay under 20 MB). The three 24 h viral systems were re-rendered with
`tools/mip_video_bleachcomp.py`: the display range follows the simulated bleaching, because with a
fixed range those movies are black after about 10 h. They are trimmed to the first 20 h (CRF 27).
