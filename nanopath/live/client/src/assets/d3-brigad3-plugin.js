// Import D3

const d3 = Object.assign({},
    require("d3")
);

// Sort functions

function sort_length(a, b) {
    return a.length - b.length;
}

function sort_index(a, b) {
    return a.index - b.index;
}


const BRIGAD3 = {
    install(Vue) {
        Vue.component('brigad3', {
            template: '<svg :height="chartHeight" :width="chartWidth" :id="chartSelector"></svg>',
            props: {
                chartData: {
                    type: Object,
                    required: true
                },  // Chart data format { ring1: [{...}, {...}] },
                chartWidth: {
                    type: Number,
                    default: 768
                },
                chartHeight : {
                    type: Number,
                    default: 768
                },
                chartSelector: {
                    type: String,
                    default: 'assemblyVisualization'  // ID of SVG element
                },
                chartScale: {
                    type: Number,
                    default: 1.0
                },
            },
            data: function () {
                return {
                    segment_color: '#F2C511',
                    segment_gap: 0.01,
                    segment_click: this.setSegmentColor,
                    segment: Object,
                    ring_gap: 20,
                    ring_width: 20,
                    inner_radius: 200,


                }
            },
            watch: {
                'chartData': {
                    handler: function () {
                        this.refreshChart();
                    },
                    deep: true
                }
            },

            computed: {
                chart_data: function () {
                    return this.chartData  // Mutable data from property
                },
                ring_sizes: function () {

                    let ring_sizes = {};

                    for (const [ring, data] of Object.entries(this.chart_data)) {
                        ring_sizes[ring] = data.reduce(function (prev, cur) {
                            return prev + cur.length;
                        }, 0);
                    }

                    return ring_sizes
                },
                ripple_size: function () {
                    let largest_ring = Object.keys(this.ring_sizes)
                        .reduce((a, b) => this.ring_sizes[a] > this.ring_sizes[b] ? a : b);

                    let ripple_size = 0;
                    let gap_size = this.segment_gap * this.ring_sizes[largest_ring];

                    this.chart_data[largest_ring].map(function (arc_segment) {
                        ripple_size += arc_segment.length + gap_size
                    });

                    return ripple_size
                },


            },
            methods: {
                initalizeChart: function () {

                    this.drawChart();
                },
                refreshChart: function () {

                    // Every-time chartData changes recompute the
                    // components ring sizes, ripple size is computed within
                    // the drawChart method

                    this.clearCanvas();
                    this.drawChart();
                },

                computeSegments: function (segments, ring_size) {

                    let i = 0;
                    let gap_size = this.segment_gap * ring_size;

                    return segments.map(function (arc_segment) {
                        arc_segment.start = i;
                        arc_segment.end = i + arc_segment.length;
                        i += arc_segment.length + gap_size
                    });


                },

                colorRings: function (color='#225ea8') {

                    // See what happens if colorRing is called
                    // is it redrawn multiple times?

                    // Color all segments of single ring in one function
                    for (const [ring, data] of Object.entries(this.chart_data)) {
                         console.log('Coloring ring: ' + ring + ' - ' + color);
                         data.map(arc_segment =>
                            arc_segment.color = color
                         );
                    }

                },

                setSegmentClick(mode) {
                    console.log('Setting segment click event.');
                    if (mode === 'recolor'){
                        this.segment_click = this.setSegmentColor
                    } else if (mode === 'information'){
                        this.segment_click = this.setSegmentColor
                    } else {
                        this.segment_click = null
                    }
                },

                setSegmentColor: function (ring, index) {

                    console.log('Set segment ' + index + ' in ring ' + ring + ' to ' + this.segment_color);

                    const segment= this.chart_data[ring].find(arc_segment =>
                            arc_segment.index === index
                    );

                    segment.color = this.segment_color;
                },


                sortRings: function (sort_mode='index') {

                    let new_data = {};

                    for (const [ring, data] of Object.entries(this.chart_data)) {
                        let ring_size = this.ring_sizes[ring];

                        if (sort_mode === 'index') {
                            // Sort by original index
                            this.computeSegments(
                                data.concat().sort(sort_index),
                                ring_size
                            )
                        } else {
                            this.computeSegments(
                                data.concat().sort(sort_length),
                                ring_size
                            )
                        }
                    }

                },

                drawChart: function () {

                    const pi = Math.PI;

                    let radiusScale = this.chartScale;
                    let segments = Object.values(this.chart_data).flat();

                    let degreeScale = d3.scaleLinear()
                        .domain([0, this.ripple_size])
                        .range([0, 360])
                    ;

                    let arc = d3.arc()
                        .innerRadius(function (d) {
                            return d.inner * radiusScale;
                        })
                        .outerRadius(function (d) {
                            return d.outer * radiusScale;
                        })
                        .startAngle(function (d) {
                            return degreeScale(d.start) * (pi / 180);
                        })
                        .endAngle(function (d) {
                            return degreeScale(d.end) * (pi / 180);
                        })
                    ;

                    console.log('Draw SVG into #' + this.chartSelector);

                    const svg = d3.select('#' + this.chartSelector);

                    const width = svg.attr('width');
                    const height = svg.attr('height');

                    let g = svg.append("g")
                        .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");

                    // let tooltip = d3.select("body").append("div")
                    //     .attr("class", "tooltip")
                    //     .style("opacity", 0);

                    // Setting up variables for possible use

                    const click_func = this.segment_click;
                    const click_color = this.segment_color;

                    let rings = g.selectAll("path")
                        .data(segments)
                        .enter()
                        .append("path")
                        .attr("d", arc)
                        .style("fill", function (d) {
                            return d.color
                        })
                        .on('click', function(d){
                            click_func(d.ring, d.index, click_color)
                        });
                    // .on('mouseover', function (d) {
                    //     tooltip.transition()
                    //         .duration(200)
                    //         .style("opacity", .9);
                    //     tooltip.html(d.text)
                    //         .style("left", (d3.event.pageX + 20) + "px")
                    //         .style("top", (d3.event.pageY + 10) + "px");
                    // })
                    // .on('mouseout', function (d) {
                    //     tooltip.transition()
                    //         .duration(200)
                    //         .style("opacity", 0)
                    // });
                },
                clearCanvas: function () {
                    d3.select("#" + this.chartSelector).selectAll("*").remove();
                }
            },
            mounted: function () {
                console.log('Mounted chart.');

                this.initalizeChart();

            },

        })
    }

};

export default BRIGAD3;
