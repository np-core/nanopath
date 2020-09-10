<template>
  <div>
    <div class="row" id="header_row">
    <div class="col-3"><gb-heading tag="h2" :style="{'color': config.secondary_color}">NanoPath <span :style="{'color': config.primary_color}">Client</span></gb-heading></div>
    <div class="col-3 header_info">
      <gb-badge class="tag" size="mini" color="turquoise">Online</gb-badge>
      <gb-badge class="tag" size="mini" color="purple">DNA</gb-badge>
      <gb-badge class="tag" size="mini" color="turquoise">SQK-RBK004</gb-badge>
    </div>
    <div class="col-5 header_info">
      <span class="run_header_title">Pipeline:</span>
      <span class="run_header_pipeline">{{ config.pipeline }}</span>
      <span class="run_header_version">{{ config.version }}</span>
      <span class="run_header_at">@</span>
      <span class="run_header_text">local, server.genomicsresearch.org</span>
      <span class="run_header_updated"> last update: </span>
      <span class="run_header_text">5 mins ago</span>
    </div>
    <div class="col-1 header_button">
      <gb-button id="select_settings" class="settings_button" left-icon="settings" circular color="turquoise" left-icon-color="#171E29" size="nano" @click="open_settings">Settings</gb-button>
    </div>
    
    <gb-divider id="header_divider" color="purple" size="large"></gb-divider>
    </div>
    <div class="row alert_row">
      <div v-if="alerts.report_group" class="alert__top"> PLEASE SELECT A SAMPLE TO GENERATE REPORT</div>
      <div v-if="alerts.report_threshold" class="alert__top"> PLEASE ENTER A POSITIVE NUMBER (READS) OR FRACTION (OF TOTAL READS) TO SET MINIMUM REPORTING THRESHOLD</div>
      <div v-if="alerts.negative_control" class="alert__top">{{ windows.report ? "EXCLUDED": "EXCLUDING" }} ALL TAXA IN NEGATIVE CONTROL ({{ selected.negative_control ? selected.negative_control.toUpperCase(): null }})</div>
    </div>

    <div class="row" id="data_row">

      <div class="row" id="main__data_window" v-show="windows.plots == true">
        <div class="col-3" id="options" >
          <div class="row" id="visual__row_select">
            <gb-select :model="selected.group" :value="selected.group" label="Sample" placeholder="Sample / Barcode" size="mini" :options="groups" class="selector" @change="select_group"></gb-select>
            <gb-select :model="selected.database" :value="selected.database" label="Database" placeholder="Database" size="mini" :options="databases" class="selector" @change="select_database" v-if="!disable_inputs"></gb-select>
            <gb-select :model="selected.negative_control" :value="selected.negative_control" label="Negative Control" placeholder="None" size="mini" :options="negative_controls" class="selector" @change="select_negative_controls"></gb-select> 
            <gb-input  :model="report.config.threshold" label="Report threshold" size="small" v-if="disable_inputs" class='selector' rounded @change="select_report_threshold"></gb-input >
            <gb-input  :model="visual.major_threshold" :value="visual.major_threshold" v-if="!disable_inputs" label="Major threshold" class="selector" rounded @change="select_major_threshold" size="small"></gb-input>
            <gb-input  :model="visual.minor_threshold" :value="visual.minor_threshold" v-if="!disable_inputs" label="Minor threshold" class="selector" rounded @change="select_minor_threshold" size="small"></gb-input>
          </div>
          <div class="row" id="visual__toggle_select">
            <gb-toggle v-model="visual.show_host" :label="settings.host_species +' (Host)'" size="small" v-if="!disable_inputs" class='toggle__toggle'></gb-toggle >
            <gb-toggle v-model="report.config.run_complete" label="Sequencing complete" size="small" v-if="disable_inputs" class='toggle__toggle'></gb-toggle >
          </div>
          <div class="row" id="visual__toggle_report">
            <gb-button class="button__report" color="black" size='medium' @click="create_report">{{ report.config.button_text }}</gb-button>
            <gb-button v-if="windows.report" class="button__return_report" left-icon="navigate_before" circular color="purple" left-icon-color="#171E29" size="nano" @click="exit_report_window"></gb-button>
          </div>
        </div>
        
        <div class="col-4" id="column__candidate_chart" v-show="!windows.report">
          <zingchart :data="major_plot.options" :series="major_plot.data"  width="480" id="plot__major" ref="plot__major" @node_click="click_major"></zingchart>
        </div>
        <div class="col-4" id="column__contamination_chart" v-show="!windows.report">
          <zingchart :data="minor_plot.options" :series="minor_plot.data"  width="480" id="plot_minor" ref="plot__minor" @node_click="click_minor"></zingchart>
        </div>
        <div class="col-8" id="window__report" v-show="windows.report">
          <div class="row">
          <div class="col-6">
          <gb-input label="Patient name" class="report__text_input" rounded v-model="report.header.patient_name"></gb-input>
          <gb-input label="Birth Date" class="report__text_input" rounded v-model="report.header.dob"></gb-input>
          <gb-input label="Location" class="report__text_input" rounded v-model="report.header.location"></gb-input>
          <gb-input label="Sample" class="report__text_input" rounded v-model="report.header.sample_type"></gb-input>
          </div>
          <div class="col-6">
          <gb-input label="Sample ID" class="report__text_input" rounded v-model="report.header.sample_id"></gb-input>
          <gb-input label="Sequenced from" class="report__text_input" rounded v-model="report.header.sequenced_from"></gb-input>
          <gb-input label="Requested by" class="report__text_input" rounded v-model="report.header.requested_by"></gb-input>
          <gb-input label="Contact" class="report__text_input" rounded v-model="report.header.contact"></gb-input>
          </div>
          </div>
        </div>
        <gb-divider id="group_divider" color="purple" size="large"></gb-divider>
          <div class="col-3">
            <div class='slide' :id="'slide_'+index" v-for="(group, index) in sorted_groups" :style="{ 'border-color': group_colors[index], 'color': group_colors[index]}" @click="select_group(group)">
              {{group}}
            </div>
          </div>
          <div class="col-9" id="window__report" v-if="windows.qc">
            <div class="row">
            <div class="col-6" id="column__plot_assembly">
            <zingchart :data="assembly_plot.options" :series="assembly_plot.data"  width="480" id="plot__assembly" ref="plot_assembly"></zingchart>
            </div>
            <div class="col-6"  id="column__plot_coverage">
            <zingchart :data="coverage_plot.options" :series="coverage_plot.data"  width="400" id="plot__coverage" ref="plot_coverage"></zingchart>
            </div>
            </div>
          </div>
        </div>
      

      <div id="main__settings_window" v-if="windows.plots == false">
        <div class="row-12">
          <span :style="{color: config.primary_color}">
          <gb-button class="button__return_settings" left-icon="navigate_before" circular color="purple" left-icon-color="#171E29" size="nano" @click="windows.plots = true">Settings</gb-button>
          Return
          </span>
        </div>
        <p class="title__settings">User interface </p>
        <div class="row" id="row__user_interface">
          <gb-toggle v-model="config.clinical_interface" label="Activate clinical interface" size="medium" id="settings__clinical_interface"></gb-toggle>
        </div>
        <p class="title__settings">Pipeline loading </p>
        <div class="row" id="row__settings">
          <gb-input label="Group regex" class="settings__setting" rounded v-model="settings.group_regex"></gb-input>
          <gb-input label="Host species" class="settings__setting" rounded v-model="settings.host_species"></gb-input>
        </div>
        
      </div>



    <gb-divider id="footer_divider" color="purple" size="large"></gb-divider>
    <div class="footer">
      <span :style="{color: config.secondary_color}">NanoPath client developed by 
      <span class="footer__author">Eike Steinig</span>, 
      <span class="footer__author">Cadhla Firth</span> and
      <span class="footer__author">Lachlan Coin</span> @
      <a href="https://github.com/np-core">https://github.com/np-core</a>
      </span>
    </div>


  </div>
  </div>
</template>



<script>

const sleep = (milliseconds) => {
  return new Promise(resolve => setTimeout(resolve, milliseconds))
}

export default {
  name: 'Dashboard',
  data() {
    return {
      config: {
        primary_color: "#AB7EF6",
        secondary_color: "#26C1C9",
        clinical_interface: false,
        pipeline: "np-pathogen",
        version: "v0.1.7"
      },
      settings: {
        group_regex: String.raw``,
        host_species : "Homo sapiens"
      },
      templates: {
        empty_data: { minor: [], major: [] },
        empty_segment: { 'background-color': null, assembly: [] },
      },
      report: {
        header: {
          patient_name: 'Anonymous',
          dob: 'not provided',
          location: 'Cairns',
          sample_type: 'Blood',
          sample_id: 'CNS-001',
          requested_by: 'Dr. Cadhla Firth',
          contact: 'cadhla.firth@jcu.edu.au',
          sequenced_from: 'Culture'
        },
        config: {
          threshold: null,
          run_complete: false,
          button_text: "Create Report"
        }
      },
      alerts: {
        negative_control: false,
        report_group: false,
        report_threshold: false
      },
      visual: {
        show_host: false,
        major_threshold: 1000,
        minor_threshold: 100
      },
      selected: {
        database: null,
        negative_control: null,
        group: null
      },
      server_data: {
        cache: {},
        reads: {},
        data: {}
      },
      windows: {
        plots: true,
        report: false,
        qc: true
      },
      disable_inputs: false,
      major_segment_selected: { 'background-color': null, assembly: [] },
      major_segments_clicked: [],
    }
  },
  computed: {

    databases: function(){
      return this.server_data.data ? [{label: "None", value: null}].concat( Object.keys(this.server_data.data).map(function(x){ return {label: x, value: x} }) ) : []
    },

    negative_controls: function(){
      if (this.server_data.data){
        if (this.selected.database){
          return [{label: "None", value: null}].concat( Object.keys(this.server_data.data[this.selected.database]).map(function(x){ return {label: x, value: x} }) )
        } else {
          return [];
        }
      } else {
        return [];
      }
    },

    groups: function(){
      if (this.server_data.data){
        if (this.selected.database){
          return [{label: "None", value: null}].concat( Object.keys(this.server_data.data[this.selected.database]).map(function(x){ return {label: x, value: x} }) )
        } else {
          return [];
        }
      } else {
        return [];
      }
    },

    selected_data: function(){     
      // Catch exception on null of server data (no results emitted from server query)
      if (this.server_data.data){
        // Catch exception on null of selected databases (no results parsed)
        if (this.selected.database){
          // Catch exception on null of selected group (no groups parsed)
          if (this.selected.group){

            var data = this.server_data.cache[this.selected.database][this.selected.group]


            console.log(this.settings.host_species)
             console.log(data)
            // Host control filter (modify selected)
            var no_human = this.filter_data(data, [this.settings.host_species], 'minor')
            console.log(no_human)
            // Negative control filter (modify selected and visualized cf. report negative control works directly of data on server)
            if (this.selected.negative_control){
              var negative_control_data = this.server_data.data[this.selected.database][this.selected.negative_control];
              var minor_negative_taxa = negative_control_data.minor.map(function(x) { return x.text });
              var major_negative_taxa = negative_control_data.major.map(function(x) { return x.text });
              var negative_control_labels = minor_negative_taxa.concat(major_negative_taxa)

              var new_minor = this.filter_data(data, negative_control_labels, 'minor');
              var new_major = this.filter_data(data, negative_control_labels, 'major');
            } else {
              var new_minor = data.minor;
              var new_major = data.major;
            }

            // Color

            var minor_n = new_minor.length;
            var major_n = new_major.length;

            if (minor_n == 1){
              var minor_scale = [this.config.secondary_color]
            } else {
              var minor_scale = this.$chroma.scale([this.config.primary_color, this.config.secondary_color]).domain([0, minor_n]).colors(minor_n);
            }
            if (major_n == 1){
              var major_scale = [this.config.primary_color]
            } else {
              var major_scale = this.$chroma.scale([this.config.primary_color, this.config.secondary_color]).domain([0, major_n]).colors(major_n);
            }
            
            new_minor.forEach(function(x,i){
                x['background-color'] = minor_scale[i]
            });
            new_major.forEach(function(x,i){
                x['background-color'] = major_scale[i]
            });

            console.log('Triggered selected data')

            console.log(new_minor, new_major)

            return { minor: new_minor, major: new_major };
          } 
        }
      }
      return this.templates.empty_data
    },
    sorted_groups: function(){
      if (this.server_data.data){
        if (this.selected.database){
          return Object.keys(this.server_data.data[this.selected.database]).sort(function (a,b) {
            return a.localeCompare(b, undefined, { numeric: true, sensitivity: 'base' });
          });
        } else {
          return [];
        }
      } else {
        return [];
      }
    },
    group_colors: function(){
      if (this.server_data.data){
        if (this.selected.database){
          const n = Object.keys(this.server_data.data[this.selected.database]).length;
          return this.$chroma.scale([this.config.primary_color, this.config.secondary_color]).domain([0, n]).colors(n);
        } else {
          return [];
        }
      } else {
        return [];
      }
    },
    server_setting: function(){
      return {
        minor_threshold: this.visual.minor_threshold,
        major_threshold: this.visual.major_threshold,
        group_regex: this.settings.group_regex,
        host_species: this.settings.host_species
      }
    },
    major_plot: function() {

      return {
        options: {
          type: 'ring',
          backgroundColor: 'transparent',
          title: { text: 'Major', color: "#26C1C9"},
          plot: {
            'border-color': '#171E29',
            'border-width': 3,
            'value-box': [
              {
                'text': "%pie-total-value\nreads",
                'placement': "center",
                'font-color': this.config.primary_color,
                'font-size': 20,
                'font-family': "Georgia",
                'font-weight': "normal"
              }
            ],
            slice: "70%",
            layout: 'auto',
            detached: false,
            tooltip: {
              text: "%t"
            }
          },
          noData: {
            text: 'Data unavailable',
            alpha: .6,
            color: "#26C1C9",
            bold: true,
            fontSize: '18px',
            textAlpha: .9
          }
        },
        data: this.selected_data.major
      }
    },
    minor_plot: function() {
      return {
        options: {
          type: 'ring',
          backgroundColor: 'transparent',
          title: { text: 'Minor', color: "#26C1C9" },
          plot: {
            'border-color': '#171E29',
            'border-width': 3,
            'value-box': [
              {
                'text': "%pie-total-value\nreads",
                'placement': "center",
                'font-color': this.config.primary_color,
                'font-size': 20,
                'font-family': "Georgia",
                'font-weight': "normal"
              }
            ],
            slice: "70%",
            layout: 'auto',
            detached: false,
            tooltip: {
                text: "%t"
              }
          },
          noData: {
            text: 'Data unavailable',
            alpha: .6,
            color: "#26C1C9",
            bold: true,
            fontSize: '18px',
            textAlpha: .9
          }
        },
        data: this.selected_data.minor
      }
    },
    coverage_plot: function() {
      return {
        options: {
          type: "area",
          backgroundColor: 'transparent',
          title: { text: 'Contig Coverage', color: "#26C1C9" },
          plot: {
            aspect: "stepped"
          },
          "scale-y": {
            guide: {
              visible: false
            }
          }
        },
        data: [
          { values: [20,40,25,50,15,45,33,34]}
        ]
      }
    },
    assembly_plot: function() {
      return {
        options: {
          type: 'ring',
          backgroundColor: 'transparent',
          title: { text: 'Metagenome Assembly', color: "#26C1C9" },
          plot: {
            'border-color': '#171E29',
            'border-width': 3,
            'value-box': [
              {
                'text': "%pie-total-value\nreads",
                'placement': "center",
                'font-color': this.major_segment_selected['background-color'],
                'font-size': 20,
                'font-family': "Georgia",
                'font-weight': "normal"
              }
            ],
            slice: "70%",
            layout: 'auto',
            detached: false,
            tooltip: {
                text: "%t"
              }
          },
          noData: {
            text: 'Data unavailable',
            alpha: .6,
            color: "#26C1C9",
            bold: true,
            fontSize: '18px',
            textAlpha: .9
          }
        },
        data: this.major_segment_selected.assembly
      }
    }
  },
  sockets: {
    get_pathogen_data(data){
      
      this.server_data.data = data.server_data.data;
      this.server_data.reads = data.server_data.reads;
    
      // Cache data unharmed as deep copy
      if (this.server_data.data){  
        this.server_data.cache = JSON.parse(
          JSON.stringify(this.server_data.data)
        );
      }

      // If there is server data (database keys)
      var databases = Object.keys(this.server_data.data);
      if (databases) {
        this.selected.database = databases[0];
      }
      if (this.sorted_groups){
        this.selected.group = this.sorted_groups[0]
      }     
    }
  },
  methods: {

    click_major: function(ev){
      console.log('Clicked major segment')
      var segment_data = this.selected_data.major[ev.plotindex]
      console.log(segment_data);
      if (!this.major_segments_clicked.includes(segment_data.text)) {
          // track clicked segments in major plot
          this.major_segments_clicked.push(segment_data.text)
      } else {
         // remove clicked segment if it is clicked again
        var _index = this.major_segments_clicked.indexOf(segment_data.text);
        console.log(_index)
        this.major_segments_clicked.splice(_index, 1);
      }

      console.log(this.major_segments_clicked)
      this.major_segment_selected = segment_data
    },
    click_minor: function(ev){
      console.log('Clicked minor', ev)
      console.log(this.selected_data.major[ev.plotindex])
    },

    // Require call to server to obtain new thresholded data (not done on client yet)
    select_major_threshold: function(threshold){
      this.visual.major_threshold = threshold;
      this.$socket.emit('get_pathogen_data', this.server_setting)    

    },
    // Require call to server to obtain new thresholded data (not done on client yet)
    select_minor_threshold: function(threshold){
      this.visual.minor_threshold = threshold;
      this.$socket.emit('get_pathogen_data', this.server_setting)    
    },
    // Triggers computed property 'selected_data' when group is changed
    select_group: function(group){
      
      this.selected.group = String.raw`${group}`;

      if (!this.selected.group){
        this.alerts.report_group = true;
      } else {
        this.alerts.report_group = false;
      }
      if (!group){
        this.report.config.threshold = null;
      }
    },
    select_report_threshold: function(user_input){

      // Handle user input as reads or percent of total reads
      
      if (isNaN(user_input)){
        // Numeric alert - input is not a number
        this.alerts.report_threshold = true
        this.report.config.threshold = null
      } else {

        if (user_input < 0){
          this.alerts.report_threshold = true
          this.report.config.threshold = null;
          return
        }

        this.alerts.report_threshold = false

        if (user_input >= 1){
          console.log('Entered number in reads', user_input)
        } else {
          console.log('Entered number in percent', user_input)
        }
        this.report.config.threshold = user_input
      }
    },
    select_database: function(db, name, e) {
      this.selected.database = db;
      console.log(this.selected_data)
    },
    select_negative_controls: function(control, name, e) {
      this.selected.negative_control = control;
      if (control){
        this.alerts.negative_control = true
      } else {
        this.alerts.negative_control = false
      }
    },
    filter_data: function(data, to_filter, plot) {
      var new_data = []
      
      console.log('Filter data 1', data)
      
      console.log('Filter data 2', data[plot])
      data[plot].map(function(x, i){
        if (!to_filter.includes(x.text)){
          new_data.push(x);
        }
      })
      console.log('Filter data 3', new_data)
      
      return(new_data)
    },
    open_settings: function() {
      this.windows.plots = false;
    },
    exit_report_window: function(){
      
      if (this.alerts.report_group){
        this.alerts.report_group = false;
      } 
      if (this.alerts.negative_control){
        this.alerts.negative_control = false;
      } 
      if (this.alerts.report_threshold){
        this.alerts.report_threshold = false;
      }

      this.windows.report = false;
      this.disable_inputs = false;
      this.selected.negative_control = null;
      if (this.sorted_groups){
        this.selected.group = this.sorted_groups[0]; // select first from sorted
      } else {
        this.selected.group = this.groups[0].value;  // always at none
      }
      
      this.report.config.threshold = null;
      this.report.config.button_text = "Create Report";
    },
    enter_report_window: function(){

      if (this.alerts.negative_control){
        this.alerts.negative_control = false;
      }
      if (this.alerts.report_group){
        this.alerts.negative_control = false;
      }

      this.windows.report = true;
      this.disable_inputs = true;
      this.selected.group = null;
      this.selected.negative_control = null;
      this.report.config.button_text = "Submit";
    },
    create_report: function(){
      if (this.windows.report){
        
        if (!this.selected.group){
          this.alerts.report_group = true
          return // cancel if no group (sample) is selected
        } else {
          this.alerts.report_group = false
        }

        if (!this.report.config.threshold){
          this.alerts.report_threshold = true
          return // cancel if no threshold is selected
        } else {
          this.alerts.report_threshold = false
        }

        if (this.selected.negative_control){
          this.alerts.negative_control = true
        } else {
          this.alerts.negative_control = false
        }

        var report_data = {};
        
        report_data['header'] = this.report.header;

        report_data['config'] = {
          pipeline: this.config.pipeline,
          version: this.config.version,
          sample: this.selected.group,
          negative_control: this.selected.negative_control,
          threshold: this.report.config.threshold,
          run_complete: this.report.config.run_complete,
          host_species: this.settings.host_species
        };

        report_data['setting'] = this.server_setting
        console.log('Report Data', report_data)
        this.$socket.emit('create_report', report_data)        
        this.exit_report_window();
      } else {
        this.enter_report_window();
      }
    },
    run_coverm: function(){
      // Method to launch CoverM against local reference genome database when selected (needs a better database selector)
      // to show quality control of genome mapping (even, uneven)

    },
    poll_server: function(){
      this.$socket.emit('get_pathogen_data', this.server_setting)
    },
    start_polling: function(){
      // Introduce slight delay to avoid potential socketio complaint
      sleep(500).then(() => {
        this.poll_server()
        setInterval(this.poll_server, 7000)
      })
    }
  },
  mounted() {

      this.poll_server();
      // this.start_polling();
  },
};
</script>


<style>

a:link {
  color: #AB7EF6;
}

a:visited {
  color: #AB7EF6;
}

#header {
  color: #26C1C9;
}
#header_client {
  color:  #AB7EF6;
}

#header_row {
  font-size: 9pt;
  margin-left: 0rem;
  margin-right: 0rem;
}

.header_info {
  margin-top: 0.5rem;
  font-size: 8pt;
}

.run_header_updated{
  margin-left: 1rem;
  color: #AB7EF6;
}

.run_header_title {
  color: #26C1C9;
}
.run_header_pipeline {
  margin-left: 1rem;
  color: #AB7EF6;
}
.run_header_version {
  margin-left: 1rem;
  color: #A9C7DF;
}
.run_header_at {
  margin-left: 1rem;
  color: #AB7EF6;
}
.run_header_text {
  margin-left: 1rem;
  color: #26C1C9;
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

#group_divider {
  margin-top: 0rem;
  margin-bottom: 2rem;
}

.footer {
  color: #26C1C9;
  font-size: 8pt;
}

#options {
  margin-left: 1rem;
}

.tag {
  margin-right: 1rem;
}

#chart__candidate {
  margin-top: 3rem;
  margin-right: 0rem;
  margin-bottom: 2rem;
}
#chart__contamination {
  margin-top: 3rem;
}

.apexcharts-datalabel-value {
  fill: #26C1C9;
}

.slide {
  padding: 0.4rem;
  margin-bottom: 0.4rem;
  width: 100%;
  font-size: 8pt;
  border-style: solid;
  border-radius: 1rem;
  background-color: #1B2431 
}
.slide:hover {
  cursor: pointer;
  background-color: #171E29;
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

.selector {
  width: 80%;
  margin-top: 1rem;
  margin-left: 1rem;
}

#select_row {
  margin-top: 0rem;
}
#visual__report_select {
  margin-top: 2rem;
}

#visual__plot_select {
  margin-top: 1rem;
}
#visual__toggle_select {
  margin-top: 1rem;
}

.button__report {
  width: 55%;
  margin-top: 1.5rem;
  margin-bottom: 1.5rem;
  margin-left: 3rem;
}

.settings_button {
  margin-top: 0.5rem;
  margin-left: 4rem;
}
.info_button {
  margin-bottom: -1rem;
  margin-left: 5rem;
}

#column__candidate_chart {
  margin-left: 0rem;
}

#column__contamination_chart {
  margin-left: 3rem;
}

.alert__top {
  padding: 0.8rem;
  margin-left: 2rem;
  margin-right: 2rem;
  margin-bottom: 0.4rem;
  width: 100%;
  font-size: 10pt;
  font-weight: bold;
  border-width: 0.1rem;
  border-style: solid;
  border-color: #26C1C9;
  color: #AB7EF6;
  border-radius: 1rem;
  text-align: center;
  background-color: #1B2431 
}

#main__settings_window {
  margin-top: 1rem;
  margin-left: 2rem;
  margin-bottom: 1rem;
  color: #26C1C9
}

.report__text_input {
  min-width: 15rem;
  margin-left: 2rem;
  margin-top: 1rem;
  margin-bottom: 1rem;
}

.settings__setting {
  min-width: 15rem;
  margin-left: 2rem;
}
#row__settings {
  margin-bottom: 2rem;
}
#settings__clinical_interface {
  margin-left: 2rem;
}
.footer__author {
  color: #AB7EF6;
}
.title__settings {
  margin-top: 2rem;
  margin-bottom: 2rem;
}
.alert_row {
  margin-left: 10rem;
}

#window__report {
  margin-top: 0.5rem;
  margin-bottom: 2rem;
}

.toggle__toggle {
  margin-top: 0.5rem;
  margin-left: 3rem;
}

.button__return_settings {
  margin-left: 1rem;
}
.button__return_report {
  margin-top: 1rem;
  margin-bottom: 1rem;
  margin-left: 10rem;
}

#column__plot_assembly{
  margin-left: 0rem;
}
#column__plot_coverage{
  margin-left: 0rem;
}
</style>