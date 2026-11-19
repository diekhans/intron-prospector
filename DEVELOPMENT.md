# Development

Notes:
* The configure script is kept current it git, which allows builds from github
  generated release download archives.

## Making a release

### Version numbering
- `AC_INIT` in `configure.ac` holds the version of the last release.
- A build from a git tree reports that version with `+main` appended, such as
  `1.5.1+main`, so a development build is never mistaken for a release.
- A build from a release archive has no `.git`, so it reports the version
  itself.
- `configure` is checked in, so GitHub release archives build without
  autoconf.

### Before tagging
- Update `CHANGELOG.md`: retitle the `next release` section to the version and
  the date, newest section first.
- Set the version in `AC_INIT` in `configure.ac`.
- `autoreconf -fi`
- `./configure`
- `make doc`, which needs pandoc.  It regenerates `src/*.man.h` from
  `docs/*.md` and the man pages in `man/`.
- `git add` the regenerated `configure`, `src/*.man.h` and `man/*.1`.

### Test
- `make distclean`
- `./configure --prefix=/some/test/dir`
- `make -j 32`
- `make -j 32 test`
- `make install`
- Run the installed `intron-prospector --version`.  A build of the git tree
  reports the release version with `+main`, since `configure` appends it
  whenever `.git` exists.
- Check the version a release archive reports, which is what users get, by
  building one that has no `.git`:
  ```
  git archive --prefix=intron-prospector-1.x.x/ v1.x.x | tar -C /some/test/dir -xf -
  cd /some/test/dir/intron-prospector-1.x.x && ./configure && make -j 32
  ./bin/intron-prospector --version
  ```
  It reports the release version with no `+main`.

### Tag and publish
- `git commit -am 'v1.x.x release'`
- `git tag v1.x.x`
- `git push`
- `git push --tags`
- On GitHub, draft a release from the tag with the CHANGELOG section for that
  version as the description.

### BioConda
- The GitHub release triggers an automatic version bump pull request against
  bioconda/bioconda-recipes.
- When its tests pass, comment on that pull request:
  `@BiocondaBot please add label`.  This clears the `Review Required` and
  `Merging is Blocked` notifications.
- Update the recipe by hand only the first time, or when the dependencies or
  the build change:
  - Fork https://github.com/bioconda/bioconda-recipes.
  - Set `version` and `sha256` in
    `recipes/intron-prospector/meta.yaml`.
  - The checksum is of the release archive:
    `curl -sL https://github.com/diekhans/intron-prospector/archive/refs/tags/v1.x.x.tar.gz | sha256sum`
  - Open a pull request against bioconda/bioconda-recipes.
- Install from the built package to check it:
  `conda create -n ip-test -c conda-forge -c bioconda intron-prospector`
