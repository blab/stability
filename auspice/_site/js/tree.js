console.log('Enter tree.js');

var freqScale = d3.scale.linear()
	.domain([0, 1])
	.range([1.5, 4.5]);

var tipRadius = 4.0;
var left_margin = 10;
var bottom_margin = 10;
var top_margin = 10;
if ((typeof branch_labels != "undefined")&&(branch_labels)) {top_margin +=15;}
var right_margin = 10;

function initDateColorDomain(intAttributes){
	var numDateValues = tips.map(function(d) {return d.num_date;})
	var minDate = d3.min(numDateValues.filter(function (d){return d!="undefined";}));
	var maxDate = d3.max(numDateValues.filter(function (d){return d!="undefined";}));	
	if (typeof time_window == "undefined"){
		time_window = maxDate-minDate;
		console.log("defining time window as " + time_window);
	} 
	if (time_window>1){
		dateColorDomain = genericDomain.map(function (d){return Math.round(10*(maxDate - (1.0-d)*time_window))/10;});
	}else{
		dateColorDomain = genericDomain.map(function (d){return Math.round(100*(maxDate - (1.0-d)*time_window))/100;});		
	}
	console.log('setting date domain '+dateColorDomain);
	dateColorScale.domain(dateColorDomain);
}


function initColorDomain(attr, tmpCS){
	//var vals = tips.filter(function(d) {return tipVisibility(d)=='visible';}).map(function(d) {return d[attr];});
	var vals = tips.map(function(d) {return d[attr];});
	var minval = d3.min(vals);
	var maxval = d3.max(vals);	
	var rangeIndex = Math.min(10, maxval - minval + 1);
	var domain = [];
	if (maxval-minval<20)
	{
		for (var i=maxval - rangeIndex + 1; i<=maxval; i+=1){domain.push(i);}
	}else{
		for (var i=1.0*minval; i<=maxval; i+=(maxval-minval)/9.0){domain.push(i);}		
	}
	tmpCS.range(colors[rangeIndex]);
	tmpCS.domain(domain);
}

function updateColorDomains(num_date){
	dateColorDomain = genericDomain.map(function(d) {return Math.round(10*(num_date - time_window*(1.0-d)))/10;});
	dateColorScale.domain(dateColorDomain);
}

function tipVisibility(d) {
	if ((d.diff < 0 || d.diff > time_window)&(date_select==true)) {
		return "hidden";
	}
	for (var k in restrictTo){
		if (d[k]!=restrictTo[k] && restrictTo[k]!="all"){
			return "hidden";
		}
	}
	return "visible";
}

function branchPoints(d) {
	var mod = 0.5 * freqScale(d.target.frequency) - freqScale(0);
	return (d.source.x-mod).toString() + "," + d.source.y.toString() + " "
		+ (d.source.x-mod).toString() + "," + d.target.y.toString() + " "
		+ (d.target.x).toString() + "," + d.target.y.toString();
}

function branchStrokeWidth(d) {
	return freqScale(d.target.frequency);
}

function branchLabelText(d) {
	var tmp_str = d.aa_muts.replace(/,/g, ', '); 
	if (tmp_str.length>50){
		return tmp_str.substring(0,45)+'...';
	}
	else {
		return tmp_str;
	}
}

function tipLabelText(d) {
	if (d.strain.length>32){
		return d.strain.substring(0,30)+'...';
	}
	else {
		return d.strain;
	}
}

function branchLabelSize(d) {
	var n = nDisplayTips;
	if (d.fullTipCount>n/15) {
		return "12px";
	}
	else {
		return "0px";
	}
}

function tipLabelSize(d) {
	if (tipVisibility(d)!="visible"){
		return 0;
	}
	var n = nDisplayTips;
	if (n<25){
		return 16;
	}else if (n<50){
		return 12;
	}else if (n<75){
		return 8;
	}
	else {
		return 0;
	}
}

function tipLabelWidth(d) {
	return tipLabelText(d).length * tipLabelSize(d) * 0.5;
}

function tree_init(){
	calcFullTipCounts(rootNode);
	calcBranchLength(rootNode);
	rootNode.branch_length= 0.01;	
	rootNode.dfreq = 0.0;
	if (typeof rootNode.pivots != "undefined"){
		time_step = rootNode.pivots[1]-rootNode.pivots[0];		
	}else{
		time_step = 1.0/12;
	}
	//setting index of frequency trajectory to use for calculating frequency change
	freq_ii = 1;
	if (typeof rootNode.pivots != "undefined") {
		if (typeof rootNode.pivots.length != "undefined") {
			freq_ii = rootNode.pivots.length - 1;
		}
	}
	calcNodeAges(LBItime_window);
	colorByTrait();
	adjust_freq_by_date();
	tree_legend = makeLegend();
	nDisplayTips = displayRoot.fullTipCount;	
}

d3.json(path + file_prefix + "tree.json", function(error, root) {

	if (error) return console.warn(error);

	nodes = tree.nodes(root);
	links = tree.links(nodes);
	var tree_legend;
	rootNode = nodes[0];
	displayRoot = rootNode;
	tips = gatherTips(rootNode, []);
	vaccines = getVaccines(tips);

	initDateColorDomain();
	if (typeof rootNode['ep'] != "undefined"){ initColorDomain('ep', epitopeColorScale);}
	if (typeof rootNode['ne'] != "undefined"){ initColorDomain('ne', nonepitopeColorScale);}
	if (typeof rootNode['rb'] != "undefined"){ initColorDomain('rb', receptorBindingColorScale);}
	date_init();
	tree_init();

	var xValues = nodes.map(function(d) {
		return +d.xvalue;
	});

	var yValues = nodes.map(function(d) {
		return +d.yvalue;
	});

	if ((typeof tip_labels != "undefined")&&(tip_labels)){
		var maxTextWidth = 0;
		var labels = treeplot.selectAll(".tipLabel")
			.data(tips)
			.enter()
			.append("text")
			.attr("class","tipLabel")
			.style("font-size", function(d) { return tipLabelSize(d)+"px"; })		
			.text(tipLabelText)
			.each(function(d) {
				var textWidth = tipLabelWidth(d);
				if (textWidth>maxTextWidth) {
					maxTextWidth = textWidth;
				}
			});
		right_margin = maxTextWidth + 10;
	}

	var xScale = d3.scale.linear()
		.domain([d3.min(xValues), d3.max(xValues)])
		.range([left_margin, treeWidth - right_margin]);

	var yScale = d3.scale.linear()
		.domain([d3.min(yValues), d3.max(yValues)])
		.range([top_margin, treeHeight - bottom_margin]);
	console.log(treeHeight +" " + top_margin);

	nodes.forEach(function (d) {
		d.x = xScale(d.xvalue);
		d.y = yScale(d.yvalue);
	});

	var clade_freq_event;
	var link = treeplot.selectAll(".link")
		.data(links)
		.enter().append("polyline")
		.attr("class", "link")
		.attr("points", branchPoints)
		.style("stroke-width", branchStrokeWidth)
		.style("stroke", branchStrokeColor)		
		.style("cursor", "pointer")
		.on('mouseover', function (d){
			linkTooltip.show(d.target, this);
			if ((colorBy!="genotype")&(typeof addClade !="undefined")){
				clade_freq_event = setTimeout(addClade, 1000, d);
			}
			})
		.on('mouseout', function(d) {
			linkTooltip.hide(d);
			if (typeof addClade !="undefined") {clearTimeout(clade_freq_event);};})		
		.on('click', function(d) {
			if ((colorBy!="genotype")&(typeof addClade !="undefined")){
				addClade(d);
			}
			displayRoot = d.target;
			nDisplayTips = displayRoot.fullTipCount;
			var dMin = 0.5 * (minimumAttribute(d.target, "xvalue", d.target.xvalue) + minimumAttribute(d.source, "xvalue", d.source.xvalue)),
				dMax = maximumAttribute(d.target, "xvalue", d.target.xvalue),
				lMin = minimumAttribute(d.target, "yvalue", d.target.yvalue),
				lMax = maximumAttribute(d.target, "yvalue", d.target.yvalue);
			if (dMax > dMin && lMax > lMin) {
				rescale(dMin, dMax, lMin, lMax);
			}
			else {
				displayRoot = d.source;
				nDisplayTips = displayRoot.fullTipCount;
				dMin = minimumAttribute(d.source, "xvalue", d.source.xvalue),
				dMax = maximumAttribute(d.source, "xvalue", d.source.xvalue),
				lMin = minimumAttribute(d.source, "yvalue", d.source.yvalue),
				lMax = maximumAttribute(d.source, "yvalue", d.source.yvalue);			
				rescale(dMin, dMax, lMin, lMax);
			}
		});

	if ((typeof branch_labels != "undefined")&&(branch_labels)){
		var mutations = treeplot.selectAll(".branchLabel")
			.data(nodes)
			.enter()
			.append("text")
			.attr("class", "branchLabel")
			.style("font-size", branchLabelSize)			
			.attr("x", function(d) {
				return d.x - 6;
			})
			.attr("y", function(d) {
				return d.y - 3;
			})
			.style("text-anchor", "end")
			.text(branchLabelText);
		}

	if ((typeof tip_labels != "undefined")&&(tip_labels)){
		treeplot.selectAll(".tipLabel").data(tips)
			.style("font-size", function(d) { return tipLabelSize(d)+"px"; })			
			.attr("x", function(d) { return d.x+10; })
			.attr("y", function(d) { return d.y+4; });			
	}

	var tipCircles = treeplot.selectAll(".tip")
		.data(tips)
		.enter()
		.append("circle")
		.attr("class", "tip")
		.attr("id", function(d) { return (d.strain).replace(/\//g, ""); })
		.attr("cx", function(d) { return d.x; })
		.attr("cy", function(d) { return d.y; })
		.attr("r", tipRadius)
		.style("visibility", tipVisibility)
		.style("fill", tipFillColor)
		.style("stroke", tipStrokeColor)
		.on('mouseover', function(d) {
			virusTooltip.show(d, this);
		})
		.on('click', function(d) {
			if ((typeof d.db != "undefined") && (d.db == "GISAID") && (typeof d.accession != "undefined")) {
				var url = "http://gisaid.org/EPI/"+d.accession;
				console.log("opening url "+url);
				var win = window.open(url, '_blank');
  				win.focus();
  			}	
  		})		
		.on('mouseout', virusTooltip.hide);

	var vaccineCircles = treeplot.selectAll(".vaccine")
		.data(vaccines)
		.enter()
		.append("text")
		.attr("class", "vaccine")
		.attr("x", function(d) {return d.x})
		.attr("y", function(d) {return d.y})
		.attr('text-anchor', 'middle')
		.attr('dominant-baseline', 'central')
		.style("font-size", "28px")
		.style('font-family', 'FontAwesome')
		.style("fill", "#555555")
		.text(function(d) { return '\uf00d'; })
		.style("cursor", "default")
		.on('mouseover', function(d) {
			virusTooltip.show(d, this);
		})
		.on('mouseout', virusTooltip.hide);


	d3.select("#reset")
		.on("click", function(d) {
			displayRoot = rootNode;
			nDisplayTips = displayRoot.fullTipCount;
			var dMin = d3.min(xValues),
				dMax = d3.max(xValues),
				lMin = d3.min(yValues),
				lMax = d3.max(yValues);
			rescale(dMin, dMax, lMin, lMax);
			removeClade();
		})

	function rescale(dMin, dMax, lMin, lMax) {

		var speed = 1500;
		
		if ((typeof tip_labels != "undefined")&&(tip_labels)){
			var maxTextWidth = 0;
			var labels = treeplot.selectAll(".tipLabel")
				.data(tips)
				.each(function(d) {
					var textWidth = tipLabelWidth(d);
					if (textWidth>maxTextWidth) {
						maxTextWidth = textWidth;
					}
				});
			right_margin = maxTextWidth + 10;
			xScale.range([left_margin, treeWidth - right_margin]);
		}		
		
		xScale.domain([dMin,dMax]);
		yScale.domain([lMin,lMax]);

		nodes.forEach(function (d) {
			d.x = xScale(d.xvalue);
			d.y = yScale(d.yvalue);
		});

		treeplot.selectAll(".tip").data(tips)
			.transition().duration(speed)
			.attr("cx", function(d) { return d.x; })
			.attr("cy", function(d) { return d.y; });

		treeplot.selectAll(".vaccine").data(vaccines)
			.transition().duration(speed)
			.attr("x", function(d) { return d.x; })
			.attr("y", function(d) { return d.y; });

		treeplot.selectAll(".link").data(links)
			.transition().duration(speed)
			.attr("points", branchPoints);

		if ((typeof tip_labels != "undefined")&&(tip_labels)){
			treeplot.selectAll(".tipLabel").data(tips)
				.transition().duration(speed)
				.style("font-size", function(d) {return tipLabelSize(d)+"px"; })			
				.attr("x", function(d) { return d.x+10; })
				.attr("y", function(d) { return d.y+4; });
		}	
			
		if ((typeof branch_labels != "undefined")&&(branch_labels)){
			console.log('shift branch_labels');
			treeplot.selectAll(".branchLabel").data(nodes)
				.transition().duration(speed)
				.style("font-size", branchLabelSize)				
				.attr("x", function(d) {  return d.x - 6;})
				.attr("y", function(d) {  return d.y - 3;});
		}

		if (typeof clades !="undefined"){
			treeplot.selectAll(".annotation").data(clades)
				.transition().duration(speed)
				.attr("x", function(d) {
					return xScale(d[1]) - 6;
				})
				.attr("y", function(d) {
					return yScale(d[2]) - 6;
				});			
		}
	}

	d3.select(window).on('resize', resize); 
	
	function resize() {
	
		containerWidth = parseInt(d3.select(".treeplot-container").style("width"), 10);
		treeWidth = containerWidth;
		treeHeight = treePlotHeight(treeWidth);
			
		d3.select("#treeplot")
			.attr("width", treeWidth)
			.attr("height", treeHeight);

		if ((typeof tip_labels != "undefined")&&(tip_labels)){
			var maxTextWidth = 0;
			var labels = treeplot.selectAll(".tipLabel")
				.data(tips)
				.each(function(d) {
					var textWidth = tipLabelWidth(d);
					if (textWidth>maxTextWidth) {
						maxTextWidth = textWidth;
					}
				});
			right_margin = maxTextWidth + 10;
		}	
			
		xScale.range([left_margin, treeWidth - right_margin]);
		yScale.range([top_margin, treeHeight - bottom_margin]);
		
		nodes.forEach(function (d) {
			d.x = xScale(d.xvalue);
			d.y = yScale(d.yvalue);
		});		
		
		treeplot.selectAll(".tip").data(tips)
			.attr("cx", function(d) { return d.x; })
			.attr("cy", function(d) { return d.y; });

		treeplot.selectAll(".vaccine").data(vaccines)
			.attr("x", function(d) { return d.x; })
			.attr("y", function(d) { return d.y; });
			
		treeplot.selectAll(".link").data(links)
			.attr("points", branchPoints);
			
		if ((typeof tip_labels != "undefined")&&(tip_labels))
		{
			treeplot.selectAll(".tipLabel").data(tips)
				.style("font-size", function(d) { return tipLabelSize(d)+"px"; })
				.attr("x", function(d) { return d.x+10; })
				.attr("y", function(d) { return d.y+4; });
		}

		if ((typeof branch_labels != "undefined")&&(branch_labels))
		{
			console.log('shift branch_labels');
			treeplot.selectAll(".branchLabel").data(nodes)
				.style("font-size", branchLabelSize)			
				.attr("x", function(d) {  return d.x - 6;})
				.attr("y", function(d) {  return d.y - 3;});
		}

		if (typeof clades !="undefined")
		{
			treeplot.selectAll(".annotation").data(clades)
				.attr("x", function(d) {
					return xScale(d[1]) - 6;
				})
				.attr("y", function(d) {
					return yScale(d[2]) - 6;
				});
		}
	}
	

	function restrictToFunc(rt) {
		restrictTo[rt] = document.getElementById(rt).value;
		console.log("restriction to "+rt+" "+restrictTo[rt]);	
		d3.selectAll(".tip")
			.style("visibility", tipVisibility);
	}

	for (rt in restrictTo){
		var tmp = document.getElementById(rt);
		if (tmp!=null){
			restrictTo[rt] = tmp.value;
		}else{restrictTo[rt]='all';}
		console.log(restrictTo);
		d3.select("#"+rt)
			.style("cursor", "pointer")
			.on("change", (function(restrictor){
							return function(){
								return restrictToFunc(restrictor);
							}
						})(rt));		
	}


	function onSelect(tip) {
		d3.select("#"+(tip.strain).replace(/\//g, ""))
			.call(function(d) {
				virusTooltip.show(tip, d[0][0]);
			});
	}

	var mc = autocomplete(document.getElementById('search'))
		.keys(tips)
		.dataField("strain")
		.placeHolder("search strains...")
		.width(800)
		.height(500)
		.onSelected(onSelect)
		.render();


	// add clade labels
	clades = rootNode["clade_annotations"];
	if (typeof clades != "undefined"){
		console.log(clades);
		var clade_annotations = treeplot.selectAll('.annotation')
			.data(clades)
			.enter()
			.append("text")
			.attr("class", "annotation")
			.attr("x", function(d) {
				return xScale(d[1]) - 6;
			})
			.attr("y", function(d) {
				return yScale(d[2]) - 6;
			})
			.style("text-anchor", "end")
			.text(function (d) {
				return d[0];
			});
		}

});

d3.json(path + file_prefix + "sequences.json", function(error, json) {
	if (error) return console.warn(error);
	cladeToSeq=json;
});

