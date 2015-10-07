function locateSequence(){
    var mutations, olap_start, olap_end;
    var seq = document.getElementById('seqinput').value.toUpperCase().replace(/[^ACGT]/g,"");
    console.log('Provided sequence:',seq);
    tmp = alignToRoot(seq);
    olap_start=tmp[0]; olap_end=tmp[1]; mutations=tmp[2];
    console.log("start, end:", olap_start, olap_end);
    console.log("mutations:", mutations);
    var bestClade = findClosestClade(mutations);
}

function findClosestClade(mutations){
    var bestClade=-1, bestScore=0;
    var tmpScore=0;
    for (clade in cladeToSeq){
        tmpScore=0;
        for (mut in mutations){
            if (stateAtPosition(clade, 'nuc', mut)==mutations[mut]){
                tmpScore++;
            }
        }
        if (clade!="root") {
            tmpScore -= 0.5*Object.keys(cladeToSeq[clade]['nuc']).length;
        }
        if (tmpScore>bestScore){
            bestScore=tmpScore;
            bestClade=clade;
        }
    }
    console.log("best match:",bestClade);
//    treeplot.selectAll('.link').filter(function(d){return d.target.clade==bestClade;})
//        .append("circle")
//        .attr("class", "userseq")
//        .attr("id","userseq")
//        .attr("r", 10)
//        .attr("cx", function(d) { console.log(d.target.strain, d.target.x); return d.target.x; })
//        .attr("cy", function(d) { return d.target.y; })
//        .style("fill", "#EE4444")
//        .style("stroke", "#EE4444");
    treeplot.selectAll('.tip').filter(function(d){return d.clade==bestClade;})
        .attr("r", function(d){console.log(d.strain); return tipRadius*2.7;})
        .style("visibility", 'visible')
        .style("fill", function (t) {
          return d3.rgb(tipFillColor(t)).brighter();
        });

}


function alignToRoot(seq){
    var rootSeq = cladeToSeq["root"]["nuc"];
    var shift = 0;
    var max_score = 0.0, max_shift;

    for(shift=0; shift<seq.length-30;shift++){
        var tmp_score = 0;
        var olaplen=Math.min(seq.length-shift, rootSeq.length);
        for (var pos=0; pos<olaplen; pos++){
            if (seq[pos+shift]==rootSeq[pos]){
                tmp_score++;
            }
        }
        tmp_score*=1.0/olaplen;
        if (tmp_score>max_score){
            max_score=tmp_score;
            max_shift=-shift;
        }
    }

    for(shift=0; shift<rootSeq.length-30;shift++){
        var tmp_score = 0;
        var olaplen=Math.min(rootSeq.length-shift, seq.length);
        for (var pos=0; pos<olaplen; pos++){
            if (seq[pos]==rootSeq[shift+pos]){
                tmp_score++;
            }
        }
        tmp_score*=1.0/olaplen;
        if (tmp_score>max_score){
            max_score=tmp_score;
            max_shift=shift;
        }
    }
    console.log("best shift: ",max_shift, " score: ",max_score);
    var mutations = {};
    if (max_shift<0){
        var olaplen=Math.min(seq.length-max_shift, rootSeq.length);
        var olap_start = 0;
        var olap_end = olaplen;
    }else{
        var olaplen=Math.min(rootSeq.length-max_shift, seq.length);
        var olap_start = max_shift;
        var olap_end = max_shift+olaplen;
    }
    for (var pos=olap_start; pos<olap_end; pos++){
        if (rootSeq[pos]!=seq[pos-max_shift]){
            mutations[pos]=seq[pos-max_shift];
        }
    }
    return [olap_start, olap_end, mutations];
}



d3.select('#seqinputsubmit').on('click', locateSequence); 
