var	vaccineChoice = {};
vaccineChoice['A/California/07/2009'] = "2009-09-25";
var vaccineStrains = Object.keys(vaccineChoice);

var restrictTo = {"region":"all"};

var genome_annotation = {'SP':[[1.2,1.2,1.2], [1,20,52]],
                         'HA1':[[1,1,1], [52,460,52+981]],
                         'HA2':[[1.2,1.2,1.2], [52+981,1200,1701]]}
