# Writing the dos

you can mix and match markdown and RST syntax in the index. I prefer markdown.

# Building the docs

To build the docs you need to have the package itself installed, you can do this by:
```bash
cd path_to_your_repo
make install
cd docs
```

Alternatively, you can install the dependencies of the package as:
```bash
cd path_to_your_repo
make python-deps
cd docs
```
and set the `PYTHONPATH` to the root folder of the repo.

To build the docs itself, do:

```bash
make build
```

To open the docs after building then do:
```bash
make open
```
