#!/usr/bin/env python

"""
Create an HTML report for the chips 1D delta tsv file. Example:

```
SORT LEFT	SORT RIGHT	GRID	NAME	NOBS	POL	k_modes	0.010409214226828742	0.011261133937174029	0.012182777180637141	0.013179850329557393	0.014258526782023065	0.015425485184587556	0.016687950783241392	0.018053740158665434	0.019531309622744485	0.021129807575989305	0.022859131150039014	0.024729987485947014	0.02675396002765626	0.02894358024112114	0.03131240520312634	0.03387510154019545	0.03664753623729868	0.039646874878603564	0.04289168792852828	0.04640206571113875	0.050199742799786756	0.05430823258715002	0.058752972868867236	0.06356148334215192	0.06876353599454123	0.07439133943774454	0.08047973832789873	0.08706642910694481	0.09419219340089238	0.10190115052006227	0.11024103062466753	0.1192634702470429	0.12902433200025415	0.13958405045256989	0.15100800630928432	0.163366931218645	0.17673734570824934	0.19120203296340235	0.2068505513808428	0.2237797890713249	0.24209456374426902	0.2619082716886795	0.28334358986850644	0.3065332354794905	0.3316207876703009	0.358761576515669	0.38812364474561284	0.41988878818531855	0.45425368134758104	0.49143109514693284	0.5316512142749634	0.5751630623933817	0.6222360439689332	0.6731616122964774	0.7282550740378158	0.787857541449114	0.8523380443841724	0.9220958151500505	0.9975627603617794	1.0792061351006974	1.1675314359335107	1.2630855307042665	1.3664600444774098	1.4782950225960478	1.5992828935353307	1.7301727560870293	1.8717750174194905	2.0249664107296237	2.1906954235538025	2.3699881703471686	2.563954745691544	2.7737960974679527	3.000811462549138	3.246406411050311	3.5121015489442757	3.7995419329233266	4.110507255800245	4.446922865511491	4.8108716859469585
0.07410502315254246	0.10306356342729851	eor0high_phase2a-128T_lst+00_19760c43	ionosub_30l_src4k_300it_8s_80kHz_i1000	20	xx	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	2.3525688E7	NaN	NaN	NaN	NaN	NaN	1.04763608E8	NaN	NaN	NaN	1.1973693E7	5.5786812E7	NaN	2074396.5	1.0956089E7	NaN	3764376.8	NaN	983839.9	8335567.0	689291.06	895753.7	304527.16	545544.8	866749.75	882131.75	542474.5	566687.25	1360905.0	2119098.0	1566835.0	2794905.2	1.5639074E7	1.38729296E8	3.03174323E9	3.13110118E9	4.49412544E8	3.7214264E7	4639440.5	3323134.5	3880235.8	8804216.0	2.08066192E8	1.87125412E10	1.61236265E10	1.27023208E8	1.4437646E7	1.8238532E7	1.42409769E10	1.36768911E10	3.112235E7	7.0417512E7	5.6278053E10	2.92188621E9	4.7905288E7	2.73219809E10	9.1176368E7	3.4383012E10	2.03877581E9	1.18260152E10	6.8236826E8	2.69318944E8	NaN	NaN	NaN	NaN
0.07410502315254246	0.10306356342729851	eor0high_phase2a-128T_lst+00_19760c43	ionosub_30l_src4k_300it_8s_80kHz_i1000	20	yy	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1.3058076E8	NaN	NaN	NaN	NaN	NaN	7.5878963E8	NaN	NaN	NaN	9.4075296E7	2.76018816E8	NaN	7468241.5	5.102868E7	NaN	1.8397906E7	NaN	2486671.5	4.5431008E7	1557785.8	1055242.2	866037.56	421840.97	615210.75	1708431.1	1684189.4	1282848.8	1839616.6	2457781.5	1627976.4	2580051.0	1.2763145E7	9.4436576E7	2.88523264E9	1.38105446E10	2.6097216E9	1.57799344E8	9205294.0	3366228.2	4420228.5	7704516.5	1.6297256E8	1.66263245E10	7.1597801E10	5.4415034E8	1.3266093E7	1.5105043E7	1.3173631E10	6.1799887E10	3.8481124E7	5.7916084E7	1.35829864E11	1.66610524E10	4.4573404E7	7.4818396E10	4.07169024E8	8.3108839E10	1.12214948E10	3.06275123E10	3.76688512E9	2.74720768E8	NaN	NaN	NaN	NaN
```

example use with singularity

```
salloc --nodes=1 --mem=20G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=6 --tmp=20G
module load singularity
mkdir /dev/shm/deleteme
cd /dev/shm/deleteme
eval singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/ssins/ssins_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/report_chips1D_tsv.py \
        --delta-tsv /astro/mwaeor/dev/nfresults/results-eor0high-ph2/chips1d_delta_lssa_fg_simple.tsv \
        --html /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/test/chips1d_delta_lssa_fg_simple_ph2.html

# in vscode host
eval singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/ssins/ssins_latest.sif python \
    -m http.server --directory /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/test/ 8000

module load singularity
singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/ssins/ssins_latest.sif bash -c "pip install d3blocks; python /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/report_chips1D_tsv.py --delta-tsv /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/test/report_chips1D.tsv --html /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/test/test.html"
```

eval singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/ssins/ssins_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/report_chips1D_tsv.py \
        --delta-tsv /astro/mwaeor/dev/nfresults/results-test/chips1d_delta_lssa_fg_simple.tsv \
        --html /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/test/chips1d_delta_lssa_fg_simple.html

"""

import pandas as pd
import argparse
import sys
import shlex
import json
import numpy as np


def get_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--delta-tsv", type=str, help="input tsv file", required=True)
    parser.add_argument("--html", type=str, help="output html file", required=True)
    return parser.parse_args(argv)

def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def render(df):
    # df = df[df['NAME'].str.contains('ionosub')].reindex()
    print(df)
    k_mode_cols = [
        col for col in
        df.columns
        if is_float(col)
    ]
    id_vars = [col for col in df.columns if col in ['GRID', 'NAME', 'NOBS', 'POL']]
    dfm = df.melt(
        id_vars=id_vars,
        value_vars=k_mode_cols,
        var_name='k_modes',
        value_name='delta'
    )
    dfm['GRID_NAME_POL'] = dfm['GRID'].astype(str) + '_' + dfm['NAME'] + '_' + dfm['POL']
    dfm['k_modes'] = dfm['k_modes'].astype(float)
    # strip commas
    dfm['delta'] = dfm['delta'].astype(str).str.replace(',', '')
    dfm['delta'] = dfm['delta'].astype(float)
    # dfm = dfm.dropna(subset=['delta'])
    print(dfm.dropna(subset=['delta']))
    # print(json.dumps(dfm.to_dict(orient='records')))
    # colFg, colFg2, colBg, colBg2 = '#657b83', '#93a1a1', '#fdf6e3', '#eee8d5'  # solarized light
    colFg, colFg2, colBg, colBg2 = '#839496', '#93a1a1', '#002b36', '#073642'  # solarized dark
    dataJson = json.dumps(dfm.to_dict(orient='records'))

    df.dropna(inplace=True, how='all', axis=1)
    for k_mode in k_mode_cols:
        if k_mode not in df:
            continue
        df[k_mode] = df[k_mode].round(0)
    df.rename(columns={col: f"{float(col):.3}" for col in k_mode_cols}, inplace=True)


    df.drop(columns=['SORT_LEFT', 'SORT_RIGHT'], inplace=True)
    df.fillna(value="", inplace=True)
    # gridsJson = json.dumps(pd.unique(dfm['GRID']).tolist())
    headerJson = json.dumps(df.columns.tolist())
    rowsJson = json.dumps(df.values.tolist())
    return f"""<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="utf-8" />
        <meta http-equiv="X-UA-Compatible" content="IE=edge" />
        <title>D3 Example</title>
        <meta name="description" content="Barchart" />
        <meta name="viewport" content="width=device-width, initial-scale=1" />
        <script src="https://cdn.polyfill.io/v2/polyfill.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.17/d3.min.js"></script>
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/css/bootstrap.min.css">
        <style>
body {{
    background-color: {colBg};
}}
        </style>
    </head>
    <body>
        <!--<h1>D3 Example</h1>-->

        <script type="module">

// adapted from:
// - https://github.com/sgratzl/d3tutorial/blob/main/examples/bare.html
// - https://observablehq.com/@d3/multi-line-chart/2
// - https://observablehq.com/@d3/line-chart-missing-data/2
// - https://d3-graph-gallery.com/graph/custom_legend.html
// - https://codepen.io/blackjacques/pen/RYVpKZ

import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";
var margin = {{ top: 40, bottom: 100, left: 100, right: 20 }};
var width = 1600 - margin.left - margin.right;
var height = 1000 - margin.top - margin.bottom;

// Creates sources <svg> element
var svg = d3
    .select("body")
    .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .attr("viewBox", [0, 0, width, height])
    .attr("style", "max-width: 100%; height: auto; overflow: visible; font: 10px sans-serif; fill: {colFg};");

// Group used to enforce margin
var g = svg.append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
svg.append("g")
    .attr("transform", "translate(" + width + margin.left + margin.right + "," + height + margin.top + margin.bottom + ")");

const data = {dataJson};

const isNumber = (x) => typeof x !== 'undefined' && x != "" && !isNaN(x);

// select valid data from the array
const validData = data.filter(d => isNumber(d.delta) && d.delta > 0);

// Create ordinal scales for grids, pols, names
const grids = validData.map(d => d.GRID).filter((v, i, a) => a.indexOf(v) === i);
const gridScale = d3.scaleOrdinal()
    .domain(grids)
    .range(d3.range(0, 1, 1/grids.length));
const gridColor = (grid) => d3.interpolateWarm(gridScale(grid));

// const pols = validData.map(d => d.POL).filter((v, i, a) => a.indexOf(v) === i);
// const polScale = d3.scaleOrdinal()
//     .domain(pols)
//     .range([0, 1]);

const namePol = (name, pol) => `${{pol}}_${{name}}`;
const namePols = validData.map(d => namePol(d.NAME, d.POL)).filter((v, i, a) => a.indexOf(v) === i);
// const dasharray = ["20,10,5,10", "20,10,5,5,5,10", "20,10,5,5,5,5,5,10", "20,10,5,5,5,5,5,5,5,10"];
const dasharray = [""];
const namePolMap = d3.scaleOrdinal()
    .domain(namePols)
    .range([0, dasharray.length-1])
    .unknown(0)
const namePolDasharray = (name, pol) => dasharray[namePolMap(namePol(name, pol))];

// Create a checkerbox selector for namePols
const namePolSelect = d3.select("body")
    .append("form")
    .attr("id", "namePolSelect")

const namePolDivs = namePolSelect.selectAll(".form-check")
    .data(namePols)
    .enter()
    .append("div")
    .attr("class", "form-check")

var namePolsSelected = namePols;

const namePolCheckboxes = namePolDivs
    .append("input")
    .attr("type", "checkbox")
    .attr("id", (d) => d)
    .attr("value", (d) => d)
    .attr("class", "form-check-input")
    .property("checked", true)
    .on("change", function(event) {{
        if (this.checked) {{
            namePolsSelected.push(this.value);
        }} else {{
            namePolsSelected = namePolsSelected.filter((d) => d !== this.value);
        }}
        console.log(namePolsSelected);
        update();
    }});

namePolDivs
    .append("label")
    .attr("for", (d) => d)
    .attr("class", "form-check-label")
    .style("color", "{colFg}")
    .text((d) => d)


// Create the positional scales.
const x = d3.scaleLog()
    .domain(d3.extent(validData, d => d.k_modes))
    .range([margin.left, width - margin.right]);

const y = d3.scaleLog()
    .domain(d3.extent(validData, d => d.delta))
    .range([height - margin.bottom, margin.top]);

// Add the horizontal axis.
svg.append("g")
    .attr("id", "xaxis")
    .attr("transform", `translate(0,${{height - margin.bottom}})`)
    .call(d3.axisBottom(x)
        .ticks(width / 80, ".2f")
    )
    .call(g => g.select(".domain").attr("stroke", "{colFg}"))
    .call(g => g.selectAll(".tick line").attr("stroke", "{colFg2}"))
    .call(g => g.selectAll(".tick text").attr("fill", "{colFg2}"))
    .call(g => g.append("text")
        .attr("x", width)
        .attr("y", margin.bottom/2)
        .attr("fill", "{colFg}")
        .attr("stroke", null)
        .attr("text-anchor", "start")
        .text("→ k(𝘩Mpc⁻¹)")
    )

// Add the vertical axis.
svg.append("g")
    .attr("id", "yaxis")
    .attr("transform", `translate(${{margin.left}},0)`)
    .attr("fill", "{colFg}")
    .call(d3.axisLeft(y)
        .ticks(height / 80, ".1e")
    )
    .call(g => g.select(".domain").attr("stroke", "{colFg}"))
    .call(g => g.selectAll(".tick text").attr("fill", "{colFg2}"))
    .call(g => g.selectAll(".tick line")
        .attr("stroke", "{colFg2}")
        .clone()
        .attr("x2", width - margin.left - margin.right)
        .attr("stroke-opacity", 0.1)
        .attr("stroke", "{colFg2}"))
    .call(g => g.append("text")
        .attr("x", -margin.left/2)
        .attr("y", 10)
        .attr("fill", "{colFg}")
        .attr("text-anchor", "start")
        .text("↑ Δ(mK²)"));

// line drawing
const line = d3.line()
    .defined(d => isNumber(d[1]))
    .curve(d3.curveStep);

const lines = svg.append("g")
    .attr("id", "lines")
    .attr("fill", "none")
    .attr("stroke-width", "1px")
    .attr("stroke-linejoin", "round")
    .attr("stroke-linecap", "round")
    .attr("stroke-opacity", 1);

// legend
// Add one dot in the legend for each name.
var size = 10
const legend = svg.append("g").attr("id", "legend");
legend.selectAll("rect")
    .data(grids)
    .enter()
    .append("rect")
        .attr("x", margin.left + width)
        .attr("y", function(d,i){{ return margin.top + i*(size+5)}}) // 100 is where the first dot appears. 25 is the distance between dots
        .style("fill", (l) => gridColor(l))
        .attr("width", size)
        .attr("height", size)
        .attr("label", (l)=>l)

// Add one dot in the legend for each name.
legend.selectAll("text")
    .data(grids)
    .enter()
    .append("text")
        .attr("x", margin.left + width + size*1.2)
        .attr("y", function(d,i){{ return margin.top + i*(size+5) + (size/2)}}) // 100 is where the first dot appears. 25 is the distance between dots
        .text( (l) => l )
        .attr("text-anchor", "left")
        .style("alignment-baseline", "middle")


// Add an invisible layer for the interactive tip.
const dot = svg.append("g")
    .attr("display", "none");

dot.append("circle")
    .attr("r", 2.5);

dot.append("text")
    .attr("text-anchor", "middle")
    .attr("y", -8)
    .attr("fill", "{colFg}");

// Add a table

const header = {headerJson};
const hdrPol = header.indexOf('POL');
if (hdrPol == -1) console.warn(`POL column not found in header: ${{header}}`);
const hdrName = header.indexOf('NAME');
if (hdrName == -1) console.warn(`NAME column not found in header: ${{header}}`);
const rows = {rowsJson};

// hashmap for each column for its color scale
const colMap = Array.from(header, (x, i) => {{
    const values = rows.map(r => r[i]).filter(isNumber);
    if (values.length === 0) return null;
    return d3.scaleLinear()
        .domain(d3.extent(values))
        .range([0, 1])
        .clamp(true);
}});

var table = d3.select("body").append("table")
    .style("border-collapse", "collapse")
    .style("border", "2px solid")
    .style("color", "{colFg}")
    .style("white-space", "nowrap");

const reFloat = new RegExp('^[0-9]+(.[0-9]+)?$');

table.append("thead").append("tr")
    .selectAll("th")
    .data(header)
    .enter().append("th")
    .text((d) => d)
    .style("text-align", (d) => reFloat.test(d) ? "right" : "left")
    .style("border", "1px {colFg} solid")
    .style("padding", "2px")
    .style("background-color", "{colBg2}")
    .style("font-weight", "bold")
    .style("font-size", "8px");

var tbody = table.append("tbody")

var points, groups, path, tr;

const update = () => {{
    // Compute the points in pixel space as [x, y, z], where z is the name of the series.
    points = data
        .filter(d => namePolsSelected.includes(namePol(d.NAME, d.POL)))
        .map((d) => [x(d.k_modes), y(d.delta), d.k_modes, d.delta, d.GRID, d.NAME, d.POL]);

    // Group the points by series.
    groups = d3.rollup(points, v => Object.assign(v, {{g: v[0][4], n: v[0][5], p: v[0][6]}}), d => `${{d[4]}}_${{d[5]}}_${{d[6]}}`);

    // console.log(Array.from(groups.keys()));

    path = lines
        .selectAll("path")
        .data(groups.values())
        .join("path")
            .attr("d", line)
            .attr("stroke", (d) => gridColor(d[0][4]))
            .attr("stroke-dasharray", (d) => namePolDasharray(d[0][5], d[0][6]))
            .attr("display", (d) => namePolsSelected.includes(namePol(d[0][5], d[0][6])) ? null : "none")

    path.exit().remove();

    var selectedRows = rows.filter((d) => namePolsSelected.includes(namePol(d[hdrName], d[hdrPol])));

    console.log(selectedRows);

    tr = tbody
        .selectAll("tr").data(selectedRows)
        .join("tr")
            .style("font-size", "8px")
            .style("padding", "2px")
    tr.exit().remove();
    tr
        .selectAll("td")
        .data((d) => d)
        .join("td")
            .style("background-color", (d, i) => colMap[i] ? d3.interpolateWarm(colMap[i](d)) : "{colBg2}")
            .text((d) => d)
            .style("text-align", (d) => reFloat.test(d) ? "right" : "left")
            .style("border", "1px {colFg} solid")

            // .attr("display", (d) => namePolsSelected.includes(namePol(d[hdrName], d[hdrPol])) ? null : "none")
}}

update();

svg
    .on("pointerenter", pointerentered)
    .on("pointermove", pointermoved)
    .on("pointerleave", pointerleft)
    .on("touchstart", event => event.preventDefault());

// When the pointer moves, find the closest point, update the interactive tip, and highlight
// the corresponding line. Note: we don't actually use Voronoi here, since an exhaustive search
// is fast enough.
function pointermoved(event) {{
    const [xm, ym] = d3.pointer(event);
    const i = d3.leastIndex(points, ([x, y]) => Math.hypot(x - xm, y - ym));
    const [x, y, k, delta, grid, name, pol] = points[i];
    pointerentered();
    path.filter(({{g}}) => g === grid)
        .style("stroke", null)
        .style("stroke-width", "2px")
        .style("stroke-opacity", null)
        .raise();
    tr.filter((d) => d[0] === grid)
        .selectAll("td")
        .style("display", null)
    dot.attr("transform", `translate(${{x}},${{y}})`);
    dot.select("text").text(`${{grid}} (${{name}}) ${{pol}}: ${{delta}} @ ${{k}}`);
    svg.property("value", data[i]).dispatch("input", {{bubbles: true}});
}}
function pointerentered() {{
    path
        .style("stroke", "{colFg2}")
        .style("stroke-opacity", 0.5)
        .style("stroke-width", "1px");
    tr
        .selectAll("td")
        .style("display", "none");
    dot.attr("display", null);
}}
function pointerleft() {{
    path.style("stroke", null)
        .style("stroke-opacity", null)
        .style("stroke-width", null);
    tr
        .selectAll("td")
        .style("display", null)
    dot.attr("display", "none");
    svg.node().value = null;
    svg.dispatch("input", {{bubbles: true}});
}}
        </script>
    </body>
</html>
"""

def main(args):
    df = pd.read_csv(args.delta_tsv, sep='\t')
    with open(args.html, 'w') as f:
        f.write(render(df))


if __name__ == '__main__':
    if len(sys.argv) > 1:
        args = get_args(sys.argv[1:])
    else:
        # is being called directly from nextflow
        args = get_args(shlex.split("""${args}"""))
    print(vars(args))
    main(args)