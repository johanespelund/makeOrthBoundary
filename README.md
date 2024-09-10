# makeOrthBoundary
OpenFOAM utility. Make cells on patches normal to the patch faces.
Tested with OpenFOAM v2206.

## Installation

```console
git clone git@github.com:johanespelund/makeOrthBoundary.git
cd makeOrthBoundary
wmake
```
This will install the binary *makeOrthBoundary* in `$FOAM_USER_APPBIN`.

## Usage
```console
makeOrthBoundary '( <patch1> ... <patchN> )' '(<excludedPatch1> ... <excludedPatchN>)'
```
