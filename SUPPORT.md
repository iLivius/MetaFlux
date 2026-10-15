# Getting help with MetaFlux

**Start with the documentation**: [iLivius.github.io/MetaFlux](https://iLivius.github.io/MetaFlux/).
The [troubleshooting page](https://iLivius.github.io/MetaFlux/troubleshooting/) describes
the failures worth knowing about in advance, each with its cause and fix; the
[configuration reference](https://iLivius.github.io/MetaFlux/reference/configuration/)
explains every key, and the per-marker pages say what to expect from each marker.

Then, depending on what you have:

| You have | Where to go |
|---|---|
| A question about how to use MetaFlux, which settings to pick, or how to read a result | [Discussions](https://github.com/iLivius/MetaFlux/discussions) |
| Something that does not work as the docs say | An issue with the **Bug report** form; [CONTRIBUTING.md](CONTRIBUTING.md) lists what to include |
| An idea for a change or a new marker | An issue with the **Feature request** form, or a discussion first if you are unsure |
| A mistake in the documentation | An issue, or a pull request if it is small |

Two things make an answer quicker: always run `snakemake -n` (a dry run) first and paste
what it prints, because the start-up checks name most config problems with their fix;
and give the MetaFlux version you run (the banner the run prints, or `git describe --tags`).

MetaFlux is maintained by one person next to a day job, so answers can take a few days.
Security-relevant problems can be reported by e-mail to the maintainer (see
[CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) for the address) rather than in a public issue.
