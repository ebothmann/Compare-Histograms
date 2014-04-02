Citing `matplotlib` in a publication
==========
See [matplotlib FAQ](http://matplotlib.org/faq/howto_faq.html#cite-matplotlib)

Configuration file format
==========

Optional top-level options
----------

- `binRange`: The first element locates the first bin (starting from 0) and
the second element gives the number of bins to be shown in the plot, like this: `[2, 12]`
- `histogramStyle`: Possible values are `bar` (default) and `step`
- `logarithmic`: Possible values are `No` (default) and `Yes`. The latter will use absolute `y` values, thus discarding the sign information.

Optional options for distributions
----------

- `isNormalized`: Possible values are `No` (default) and `Yes`. If `Yes`, all bin heights will be normalized to the ones in this distribution. Only one distribution may have this option set to `Yes`

Optional options for a YODA distribution
----------

- `histogramIndex`: YODA files may consist of several histograms. This option specifies which one is to be used. If the option is not given, the first histogram with enough bins present will be used

Optional options for a non-YODA distribution
----------

- `column`: Gives the values column in the text file to be used as bin heights. Column numbering starts at 0. If this option is not given, the last column will be used