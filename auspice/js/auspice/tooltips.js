var virusTooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0, 12])
	.html(function(d) {
	
		string = "";
				
		// safe to assume the following attributes
		if (typeof d.strain != "undefined") {
			string += d.strain;
		}
		string += "<div class=\"smallspacer\"></div>";
		
		string += "<div class=\"smallnote\">";		
		
		// check if vaccine strain
		if (vaccineStrains.indexOf(d.strain) != -1) {
			string += "Vaccine strain<br>";
			var vaccine_date = new Date(vaccineChoice[d.strain]);

			string += "First chosen " + vaccine_date.toLocaleString("en-us", { month: "short" }) + " " + vaccine_date.getFullYear() + "<br>";
			string += "<div class=\"smallspacer\"></div>";
		}			
		
		if (typeof d.country != "undefined") {
			string += d.country.replace(/([A-Z])/g, ' $1');
		}
		if (typeof d.date != "undefined") {
			string += ", " + d.date;
		}
		if ((typeof d.db != "undefined") && (typeof d.accession != "undefined") && (d.db == "GISAID")) {
			string += "<br>GISAID ID: EPI" + d.accession;
		}
		if (typeof d.lab != "undefined") {
			if (d.lab != "") {
				string += "<br>Source: " + d.lab.substring(0,25);
				if (d.lab.length>25) string += '...';
			}
		}			
		string += "</div>";
		
		string += "<div class=\"smallspacer\"></div>";
				
		// following may or may not be present
		string += "<div class=\"smallnote\">";
		if (typeof d.ep != "undefined") {
			string += "Epitope distance: " + d.ep + "<br>";
		}
		if (typeof d.ne != "undefined") {
			string += "Non-epitope distance: " + d.ne + "<br>";
		}
		if (typeof d.rb != "undefined") {
			string += "Receptor binding distance: " + d.rb + "<br>";
		}
		if (typeof d.LBI != "undefined") {
			string += "Local branching index: " + d.LBI.toFixed(3) + "<br>";
		}
		string += "</div>";
		return string;
	});
treeplot.call(virusTooltip);


var linkTooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0, 12])
	.html(function(d) {
		string = ""
		if (typeof d.frequency != "undefined") {
			string += "Frequency: " + (100 * d.frequency).toFixed(1) + "%"
		}
		string += "<div class=\"smallspacer\"></div>";
		string += "<div class=\"smallnote\">";
		if ((typeof d.aa_muts !="undefined")&&(mutType=='aa')){
			var ncount = 0;
			for (tmp_gene in d.aa_muts) {ncount+=d.aa_muts[tmp_gene].length;}
			if (ncount) {string += "<b>Mutations:</b><ul>";}
			for (tmp_gene in d.aa_muts){
				if (d.aa_muts[tmp_gene].length){
					string+="<li>"+tmp_gene+":</b> "+d.aa_muts[tmp_gene].replace(/,/g, ', ') + "</li>";
				}
			}
		}
		else if ((typeof d.nuc_muts !="undefined")&&(mutType=='nuc')&&(d.nuc_muts.length)){
			var tmp_muts = d.nuc_muts.split(',');
			var nmuts = tmp_muts.length;
			tmp_muts = tmp_muts.slice(0,Math.min(10, nmuts))
			string += "<li>"+tmp_muts.join(', ');
			if (nmuts>10) {string+=' + '+ (nmuts-10) + ' more';}
			string += "</li>";
		}
		string += "</ul>";
		string += "click to zoom into clade"
		string += "</div>";
		return string;
	});
treeplot.call(linkTooltip);


var matchTooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0, 12])
	.html(function(d) {
		string = d.strain+ "<i> is closest match of:</i><ul>";
		string += "<div class=\"smallspacer\"></div>";
		for (var mi=0; mi<d.matches.length;mi++){
			string+="<li>" +d.matches[mi].substring(0,Math.min(30,d.matches[mi].length))+'</li>';
		}
		string += "</ul>";
		return string;
	});
treeplot.call(matchTooltip);
