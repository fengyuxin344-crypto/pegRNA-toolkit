# Starting the toolkit

Every time you use it:

```bash
conda activate snake_libanalysis
cd ~/Documents/pegRNA/pegRNA-toolkit
caffeinate -dimsu python app.py
```

It prints `running at http://localhost:5050` and opens your browser. If it
does not open, go to <http://localhost:5050> manually.

Stop it with **Control + C**.

!!! warning "Activate `snake_libanalysis` first"
    If you start the app from a different environment, the Analysis tab
    cannot find snakemake.

## Why `caffeinate`

It keeps the Mac awake through long runs. Both PRIDICT prediction and NGS
analysis can take a couple of hours, and a sleeping Mac interrupts them.

## Updating the toolkit

You do not need a new copy of the folder when something changes:

```bash
cd ~/Documents/pegRNA/pegRNA-toolkit
git pull
```

That is the whole update. Your own files are untouched — `optiprime-src/`,
run outputs and local data are excluded from the repository, so `git pull`
only ever changes the toolkit code.

If `git pull` complains about local changes, you have edited the toolkit
itself. Either keep your version:

```bash
git stash
git pull
git stash pop
```

or discard it:

```bash
git checkout -- <file>
git pull
```

## Terminal basics

| | |
|---|---|
| `cd` | change directory |
| **Control + C** | stop the running command |
| **↑** | recall the previous command |
| **Cmd + K** | clear the screen |
