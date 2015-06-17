## Analysis of stabilizing and destabilizing mutations in influenza evolution

## Pipeline

### Generate list of mutations

From `augur/` directory run nextflu pipeline to filter, align, build a tree and do ancestral state reconstruction on this tree:

```
python src/H3N2_process.py -v 1 -y 10 --stop refine
```

Run `tree_refine` to annotate tree with trunk vs side branch etc and list mutations:

```
python src/H3N2_process.py --start mutations --stop mutations
```



## License and copyright

This software is based on the nextflu pipeline.

nextflu is copyright 2014-2015 Trevor Bedford and Richard Neher.

Source code to nextflu is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). nextflu is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
