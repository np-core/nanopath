<template>
  <div>
  <div class="row" id="header_row">
      <div class="col-3"><gb-heading tag="h2">NanoPath - Alive - v0.1.0</gb-heading></div>>
      
      <div class="col-3 header_info">
      <gb-badge class="tag" size="mini" color="turquoise">Online</gb-badge>
      <gb-badge class="tag" size="mini" color="purple">DNA</gb-badge>
      <gb-badge class="tag" size="mini" color="turquoise">SQK-RPB004</gb-badge>
      </div>
      <div class="col-4 header_info"><span class="run_header_title">Pipeline: </span><span class="run_header_text">{{ config.pipeline }} - {{ config.run_id }} - {{ config.date }}</span></div>
      <gb-toggle  id="serverlog_toggle" class="header_toggle" size="mini" label="Serverlog" v-model="toggle_serverlog"></gb-toggle>
      <gb-toggle  id="config_toggle" class="header_toggle" size="mini" label="Configs" v-model="toggle_config"></gb-toggle>
      <gb-divider id="header_divider" color="purple" size="large"></gb-divider>

  </div>
  <div class="row" id="data_row">

    <div class="col-3" id="data_display">

      <div id="overview">
      
      <div id="toggles">
        <div class="row">
          <div class="col-4" id="toggles_left">
          <gb-toggle id="barcode_toggle" class="dashboard_toggle_1" size="mini" label="Barcodes" v-model="toggle_barcodes"></gb-toggle>
          </div>
          <div class="col-5" id="toggles_middle">
          <gb-toggle id="analysis_toggle" class="dashboard_toggle_1" size="mini" label="Analysis" v-model="toggle_analysis"></gb-toggle>
          </div>
        </div>
        <div class="row">
        <div class="col-4" id="toggles_left">
          <gb-toggle id="achilles_toggle" class="dashboard_toggle_2" size="mini" label="Achilles" v-model="toggle_achilles"></gb-toggle>
        </div>
        <div class="col-4" id="toggles_middle">
          <gb-toggle id="sketchy_toggle" class="dashboard_toggle_2" size="mini" label="Sketchy" v-model="toggle_sketchy"></gb-toggle>
        </div>
        </div>
      </div>

      <gb-divider id="body_divider" color="purple" size="large"></gb-divider>
      
      <div id="summary">
        <p id="container"><span class="run_header_title">Basecalled: </span><span class="summary_live">{{ online.basecalled }}</span><p>
        <p id="container"><span class="run_header_title">Reads: </span><span class="summary_live">{{ new Intl.NumberFormat().format(online.reads) }}</span></p>
        <p id="container"><span class="run_header_title">Basepairs: </span><span class="summary_live">{{ new Intl.NumberFormat().format(online.bp) }}</span></p>
      </div>

      
      <div id="config">  
        <p id="container"><span class="run_header_title">Read N50: </span><span class="summary_config">{{ new Intl.NumberFormat().format(online.n50) }} bp</span></p> 
        <p id="container"><span class="run_header_title">Mean length: </span><span class="summary_config">{{ new Intl.NumberFormat().format(online.mean_length) }} bp</span></p>    
        <p id="container"><span class="run_header_title">Median length: </span><span class="summary_config">{{ new Intl.NumberFormat().format(online.median_length) }} bp</span></p>    
        <p id="container"><span class="run_header_title">Mean quality: </span><span class="summary_config">{{ online.mean_quality }}</span></p>    
        <p id="container"><span class="run_header_title">Median quality: </span><span class="summary_config">{{ online.median_quality }}</span></p>   
      </div>

      <gb-divider id="body_divider" color="purple" size="large"></gb-divider>
      
      <p id="container"><span class="run_header_title">Docker: </span><span id="container_text">{{ config.container }} @ {{ config.host }}</span></p>
      <p id="container"><span class="run_header_title">Device: </span><span id="container_text">{{ config.device_id }}, {{ config.device_type }}</span></p>
      <p id="container"><span class="run_header_title">Flowcell: </span><span id="container_text">{{ config.flow_cell_id }}, {{ config.flow_cell_product_code }}</span></p>
      <p id="container"><span class="run_header_title">Run: </span><span id="container_text">{{ config.run_id }}, {{ config.sample_id }}</span></p>

      
      <gb-divider id="footer_divider" color="purple" size="large"></gb-divider>
      <span id="display_footer"><a href="https://github.com/np-core">https://github.com/np-core</a></span>

      </div>
      
     </div>
    
    <div class="col-9" id="data_dashboard">
      <div id="serverlog_timeline" v-if="toggle_serverlog">
        <vue-timeline-update
            :date="new Date('2016-09-30')"
            title="Pipeline"
            description="Test: np-sepsis v0.3.1"
            category="pipeline"
            icon="code"
            color="purple"
            is-last
          />
          <vue-timeline-update
            :date="new Date('2016-09-30')"
            title="Pipeline"
            description="Started: np-sepsis v0.3.1"
            category="pipeline"
            icon="code"
            color="purple"
            size="small"
            is-last
          />
      </div>
      <div id="run_charts" v-if="!toggle_serverlog">
        <gb-tabs class="barcode_tabs" :tabs="tabs" full-width="true" size="mini" v-if="toggle_barcodes"></gb-tabs>
        <div class="row" id="run_row_1">
        <div class="col-6">
        <apexchart width="450" height="300" type="bar" :options="read_length_chart.options" :series="read_length_chart.series"></apexchart>
        </div>
        <div class="col-6">
        <apexchart width="450" height="300" type="bar" :options="read_quality_chart.options" :series="read_quality_chart.series"></apexchart>
        </div>
        </div>
        <div class="row" id="run_row_2">
        <div class="col-6">
        <apexchart width="450" height="300" type="bar" :options="barcode_chart.options" :series="barcode_chart.series"></apexchart>
        </div>
        <div class="col-6">
        <apexchart width="450" height="300" type="area" :options="yield_area_chart.options" :series="chart_data.yield_series"></apexchart>
        </div>
        <div class="col-6">
        
        </div>
        </div>
      </div>
    </div>
  </div>

  </div>
</template>

<style>

a:link {
  color: #A9C7DF;
}

a:visited {
  color: #A9C7DF;
}

#header {
  color: #A9C7DF;
}

#overview {
  color: #A9C7DF;
}

#run_charts {
  color: #A9C7DF;
}

.run_header_title {
  color: #A9C7DF;
}
.run_header_text {
  color: #26C1C9;
}

#header_row {
  font-size: 9pt;
  margin-left: 0rem;
  margin-right: 0rem;
}

.header_info {
  margin-top: 0.5rem;
}

#header_divider {
  margin-top: 1rem;
  margin-bottom: 1rem;
}

#body_divider {
  margin-top: 0.5rem;
  margin-bottom: 0.5rem;
}

#footer_divider {
  margin-top: 1rem;
  margin-bottom: 0.1rem;
}

#run_id {
  margin-top: 2rem;
  margin-bottom: 1rem;
  color: #26C1C9
}

#container  {
  margin-top: 1rem;
  font-size: 10pt;
  color: #AB7EF6;
}

#container_text {
  float: right;
}

#summary {
  margin-top: 1rem;
}

#toggles {
  margin-bottom: 1rem;
  margin-top: 1rem;
}

.header_toggle {
  margin-right: 1.5rem;
  margin-top: 0.7rem;
  zoom: 0.8;
  -moz-transform: scale(0.8);
}

#display_footer {
  font-size: 8pt;
}

.tag {
  margin-right: 1rem;
}

.dashboard_toggle_1 {
  margin-bottom: 0.8rem;
}

.summary_live {
    color: #26C1C9;
    float: right;
    font-size: 10pt;
}

.summary_config {
    color: #A3E9EF;
    float: right;
    font-size: 10pt;
}

.barcode_tabs {
  margin-bottom: 2rem;
}



</style>


<script>

const range = (start, stop, step = 1) =>
  Array(Math.ceil((stop - start) / step)).fill(start).map((x, y) => x + y * step)

export default {
  name: 'Dashboard',
  data() {
    return {

        toggle_barcodes: false,
        toggle_analysis: false,
        toggle_sketchy: false,
        toggle_achilles: false,
        toggle_serverlog: false,
        toggle_config: false,

        online: {
          basecalled: '100 %',
          reads: 1092738,
          bp: 87898372,
          n50: 17876,
          mean_length: 5872,
          median_length: 7625,
          mean_quality: 10.12,
          median_quality:13.76
        },

        barcodes: [
          { label: '1', reads: 6253, test: [ {x: new Date('2020-04-03T15:33Z').getTime(), y: 0}, {x: new Date('2020-04-03T16:44Z').getTime(), y: 50}] },
          { label: '2', reads: 3253, test: [ {x: new Date('2020-04-03T15:33Z').getTime(), y: 0}, {x: new Date('2020-04-03T16:44Z').getTime(), y: 80}] },
          { label: '3', reads: 3253, test: [ {x: new Date('2020-04-03T15:33Z').getTime(), y: 0}, {x: new Date('2020-04-03T16:44Z').getTime(), y: 80}] },
          { label: '4', reads: 3253, test: [ {x: new Date('2020-04-03T15:33Z').getTime(), y: 0}, {x: new Date('2020-04-03T16:44Z').getTime(), y: 80}] },
          { label: '5', reads: 3253, test: [ {x: new Date('2020-04-03T15:33Z').getTime(), y: 0}, {x: new Date('2020-04-03T16:44Z').getTime(), y: 80}] },
          { label: '6', reads: 3253, test: [ {x: new Date('2020-04-03T15:33Z').getTime(), y: 0}, {x: new Date('2020-04-03T16:44Z').getTime(), y: 80}] },
          { label: '7', reads: 3253, test: [ {x: new Date('2020-04-03T15:33Z').getTime(), y: 0}, {x: new Date('2020-04-03T16:44Z').getTime(), y: 80}] },
          { label: '8', reads: 6253, test: [ {x: new Date('2020-04-03T15:33Z').getTime(), y: 0}, {x: new Date('2020-04-03T16:44Z').getTime(), y: 50}] },
          { label: '9', reads: 3253, test: [ {x: new Date('2020-04-03T15:33Z').getTime(), y: 0}, {x: new Date('2020-04-03T16:44Z').getTime(), y: 80}] },
          { label: '10', reads: 3253, test: [ {x: new Date('2020-04-03T15:33Z').getTime(), y: 0}, {x: new Date('2020-04-03T16:44Z').getTime(), y: 80}] },
          { label: '11', reads: 3253, test: [ {x: new Date('2020-04-03T15:33Z').getTime(), y: 0}, {x: new Date('2020-04-03T16:44Z').getTime(), y: 80}] },
          { label: '12', reads: 3253, test: [ {x: new Date('2020-04-03T15:33Z').getTime(), y: 0}, {x: new Date('2020-04-03T16:44Z').getTime(), y: 80}] }

        ],
        
        config: {
          directory: '/data/ont/runs/7c254009/fast5_pass',
          sample_id: '20200323_Saureus_MP',
          run_id: '7c254009',
          host: 'aithm-server',
          gpu: 'NVIDIA GTX-1080 Ti',
          uuid: '7c254009-422b-46a2-8c69-d8a00d97363b',
          container: 'np-live:latest',
          pipeline: 'np-sepsis v0.3.1',
          date: '20-04-2020 16:20:00',
          basecaller: 'Guppy v3.3.0',
          model: 'dna_r9.4.1_hac_450bps',
          flow_cell_id: 'FAK11919',
          flow_cell_product_code: 'FLO-FLG001',
          device_id: 'MN17546',
          device_type: 'MINION',
          min_length: 200,
          min_quality: 7,
        },

        chart_settings: {
          barcode_line: 5000,
          read_length_bins: range(5000, 85000, 5000),
          read_quality_bins: range(0, 21, 1)
        },

        chart_data: {
          read_length_distribution: [27873, 24873, 23873, 23873, 20873, 19873, 18873, 17873, 14873, 12873, 10873, 7873, 6873, 5873, 4873, 2873],
          read_quality_distribution: [0, 0, 0, 0, 0, 3, 4, 76, 128, 245, 432, 765, 543, 426, 253, 176, 86, 22, 0, 0, 0],
          yield_series: []
        },

        user_settings: {
          read_distributions: {
            length_bin_min: 0,
            length_bin_max: 80000,
            length_bin_size: 1000,
            quality_bin_min: 0.0,
            quality_bin_max: 21.0,
            quality_bin_size: 1.0
          }
        }

        

    }
  },
  computed: {

    tabs: function() {
      return this.barcodes.map(x => ({label: x.label, value: x.label}));
    },

    barcode_chart: function() {
      return {
        options: {
          chart: {
            id: 'barcode-read-chart',
            toolbar: {
              show: false
            }
          },
          title: {
            text: "Barcode reads",
            align: 'center',
          },
          annotations: {
            yaxis: [ { y: this.chart_settings.barcode_line } ]
          },
          grid: {
            show: false
          },
          legend: {
            show: false
          },
          tooltip: {
            enabled: true,
            theme: "dark",
            fillSeriesColor: true
          },
          plotOptions: {
            bar: {
              distributed: true
            }
          },
          dataLabels: {
            enabled: false
          },
          xaxis: {
            categories: this.barcodes.map(x => x.label),
            labels: {
              style: {
                  fontSize: '10px'
              }
            }
          },
          yaxis: {
            labels: {
              style: {
                  fontSize: '10px'
              }
            }
          },
          colors: this.$chroma.scale(['#AB7EF6', '#26C1C9']).domain([0, this.barcodes.length]).colors(this.barcodes.length),

        },
        series: [{
          name: 'Reads',
          data: this.barcodes.map(x => x.reads)
        }]
      }
    },

    read_length_chart: function() {
      return {
        options: {
          chart: {
            id: 'read-length-chart',
            toolbar: {
              show: false
            }
          },
          title: {
            text: "Read length (bp)",
            align: 'center',
          },
          grid: {
            show: false
          },
          legend: {
            show: false
          },
           tooltip: {
            enabled: true,
            theme: "dark",
            fillSeriesColor: true
          },
          plotOptions: {
            bar: {
              distributed: true
            }
          },
          dataLabels: {
            enabled: false
          },
          xaxis: {
            categories: this.chart_settings.read_length_bins,
            labels: {
              style: {
                  fontSize: '10px'
              }
            }
          },
          yaxis: {
            labels: {
              style: {
                  fontSize: '10px'
              }
            }
          },
          colors: this.$chroma.scale(['#AB7EF6', '#26C1C9']).domain([0, this.chart_settings.read_length_bins.length]).colors(this.chart_settings.read_length_bins.length),

        },
        series: [{
          name: 'Reads',
          data: this.chart_data.read_length_distribution
        }]
      }
    },

    read_quality_chart: function() {
      return {
        options: {
          chart: {
            id: 'read-quality-chart',
            toolbar: {
              show: false
            }
          },
          title: {
            text: "Read quality (Q)",
            align: 'center',
          },
          grid: {
            show: false
          },
          legend: {
            show: false
          },
           tooltip: {
            enabled: true,
            theme: "dark",
            fillSeriesColor: true
          },
          plotOptions: {
            bar: {
              distributed: true
            }
          },
          dataLabels: {
            enabled: false
          },
          xaxis: {
            categories: this.chart_settings.read_quality_bins,
            labels: {
              style: {
                  fontSize: '10px'
              }
            }
          },
          yaxis: {
            labels: {
              style: {
                  fontSize: '10px'
              }
            }
          },
          colors: this.$chroma.scale(['#AB7EF6', '#26C1C9']).domain([0, this.chart_settings.read_quality_bins.length]).colors(this.chart_settings.read_quality_bins.length),

        },
        series: [{
          name: 'Reads',
          data: this.chart_data.read_quality_distribution
        }]
      }
    },

    yield_area_chart: function(){
      return {
        options: {
          chart: {
            type: 'area',
            stacked: true,
            toolbar: {
              show: false
            }
          },
           title: {
            text: "Barcode yields (Mbp)",
            align: 'center',
          },
          colors: this.$chroma.scale(['#AB7EF6', '#26C1C9']).domain([0, this.barcodes.length]).colors(this.barcodes.length),
          grid: {
              show: false
          },
          legend: {
            show: false
          },
          fill: {
            gradient: {
              shade: 'dark',
              shadeIntensity: 0.8,
            }
          },
          dataLabels: {
            enabled: false
          },
          stroke: {
            curve: 'smooth',
            width: 1
          },
          xaxis: {
            type: 'datetime',
            labels: {
              style: {
                  fontSize: '10px'
              }
            }
          },
          yaxis: {
            labels: {
              style: {
                  fontSize: '10px'
              }
            }
          },
        },
      }
    },
   
  },
  sockets: {
    get_live_update(data){

      console.log(data);
    }

  },
  methods: {
    get_barcode_yields: function() {
      return this.barcodes.map(function(x){
        return new Object({name: x.label, data: x.test})
      });
    }
  },
  mounted() {
    // For some reason the yield area chart crashes with a computed property assigned to series
    // This is a work-around for now since any update to the data will be triggered by socker transfer anyway.
    this.$socket.emit('get_live_update', {user_settings: this.user_settings})
    this.chart_data.yield_series = this.get_barcode_yields();  
  },
};
</script>

