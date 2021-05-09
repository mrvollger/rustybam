var margin = {top: 20, right: 20, bottom: 40, left: 40};

d3.tsv("data.tsv", function(data) {
    for (var i = 0; i < data.length; i++) {
        console.log(data[i].query_name);
    }
});

var aln_data = [
          {c1_nm: "c1", c1_st: 0, c1_en: 100, c2_nm: "c2", c2_st: 100, c2_en: 200},
          {c1_nm: "c1", c1_st: 100, c1_en: 300, c2_nm: "c2", c2_st: 400, c2_en: 200},
          {c1_nm: "c1", c1_st: 300, c1_en: 400, c2_nm: "c2", c2_st: 100, c2_en: 200},
          //{c1_nm: "c1", c1_st: 300, c1_en: 400, c2_nm: "c3", c2_st: 300, c2_en: 200}
      ];
//var ct_names = d3.set(aln_data, function(d){return d.c2_nm;});
var t_names = [...new Set(aln_data.map(d => d.c1_nm))];
var q_names = [...new Set(aln_data.map(d => d.c2_nm))];
var ct_names = q_names.concat(t_names);
console.log(ct_names);

var height=100 * (ct_names.length);//*d3.set(ct_names).size();
var width=800;
   
 var xscale = d3.scaleLinear()
          .domain([d3.min(aln_data, function(d) { return d3.min([d.c1_st,d.c2_st]) }),
                   d3.max(aln_data, function(d) { return d3.max([d.c1_en,d.c2_en]) })])
          .range([margin.left, width - margin.right])
          .nice();

 var yscale_d = d3.scaleBand()
          .domain(ct_names)
          .range([height - margin.bottom, margin.top])
          .paddingInner(1)
          .align(0);

var container = d3.select("#chart")
    .append("svg")
    .attr("width", "100%")
    //.attr("height", "100%")
    .attr("viewBox", `0 0 ${width} ${height}`); // top, left, width, down




// my function 
function help_draw_alignment(c1_nm, o_c1_st, o_c1_en, c2_nm, o_c2_st, o_c2_en) {
    
    var c1_st = xscale(o_c1_st), c1_en = xscale(o_c1_en),
    c2_st = xscale(o_c2_st), c2_en = xscale(o_c2_en);
    
    const path = d3.path(),
    c1_h = yscale_d(c1_nm),//+yscale_d.bandwidth(),
    c2_h = yscale_d(c2_nm), //yscale_d(c2_nm),
    mid = (c1_h + c2_h) / 2; //yscale((c1_h+c2_h)/2);
    
    // connect c1 start and end
    path.moveTo(c1_st, c1_h);
    path.lineTo(c1_en, c1_h);
    // connect the ends ends
    path.bezierCurveTo(c1_en, mid, c2_en, mid, c2_en, c2_h);
    // at contig 2 end go to c2 start 
    path.lineTo(c2_st, c2_h);
    // make a bezier the goes from c2_st to c1_st 
    path.bezierCurveTo(c2_st, mid, c1_st, mid, c1_st, c1_h);
    // path.bezierCurveTo(cpx1, cpy1, cpx2, cpy2, x1, y1);
    path.closePath();

    var color = "#af0404" // red
    if( c2_en < c2_st){
    var color = "#3282b8";
    }
    
    container.append("line")
        .attr("x1",c1_st+10).attr("y1",c1_h)
        .attr("x2",c1_en-5).attr("y2",c1_h)
        .attr("stroke", "black")  
        //.attr("stroke-width",2)
        .attr("marker-start","url(#arrow)")
        .attr("marker-end","url(#arrow)");
    // target text         
    container.append('text')
        .attr("x",(c1_st+c1_en)/2).attr("y",c1_h+10)
        .style("fill", "black")
        .style("font-size", "10px")
        .attr("text-anchor", "middle")
        .attr("font-weight", "normal") 
        .text(`${o_c1_st} - ${o_c1_en}`);
    // query text
    container.append('text')
        .attr("x",(c2_st+c2_en)/2).attr("y",c2_h-10)
        .style("fill", "black")
        .style("font-size", "10px")
        .attr("text-anchor", "middle")
        .attr("font-weight", "normal") 
        .text(`${o_c2_st} - ${o_c2_en}`);

    // make the highlight regions 
    container.append("path")
        .attr("d", path)
        .attr("stroke", "black")
        .attr("fill", color)
        .attr('opacity', '.4')
        .on('mouseover', function () {
            d3.select(this).transition()
                    .duration(1)
                    .attr('opacity', '.65');
        })
        .on('mouseout', function () {
            d3.select(this).transition()
                    .duration(1)
                    .attr('opacity', '.4');
        })
    
}
// format the d as input for drawing the alignment
function draw_alignment(d, i){
    help_draw_alignment(d.c1_nm, d.c1_st, d.c1_en, d.c2_nm, d.c2_st, d.c2_en);
}

// draw the x axis 
container.append('g')
    .attr('transform', `translate(0, ${10+height-margin.bottom})`)
    .call(d3.axisBottom(xscale));

// draw the y axis 
container.append('g')
    .attr('transform', `translate(${margin.left}, 0)`)
    .call(d3.axisLeft(yscale_d));

// add in the data 
container.selectAll('g.item')
    .data(aln_data)
    .enter()
    .each(draw_alignment)
    .selectAll('path');
