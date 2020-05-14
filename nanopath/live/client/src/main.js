import Vue from 'vue';
import VueDarkMode from "@growthbunker/vuedarkmode";

import App from './App';
import router from './router';
import BootstrapVue from 'bootstrap-vue';

Vue.use(VueDarkMode, {
    theme: "dark",
    component: [
        "alert", "avatar", "badge", "button", "divider", "heading", "icon",  "progress-bar",  "spinner",
        "checkbox", "file", "input", "label", "message", "radios", "select", "tabs", "textarea", "toggle"
    ]
});
Vue.use(BootstrapVue);

import 'bootstrap/dist/css/bootstrap.css'
import 'bootstrap-vue/dist/bootstrap-vue.css'

// Data Viz

import chroma from 'chroma-js';

Object.defineProperty(Vue.prototype, '$chroma', { value: chroma });

import TaxDonut from './assets/d3-taxcomp-plugin';

Vue.use(TaxDonut);

import VueApexCharts from 'vue-apexcharts'
Vue.use(VueApexCharts);

Vue.component('apexchart', VueApexCharts);

import vuetimeline from "@growthbunker/vuetimeline"

Vue.use(vuetimeline);

// Fontawesome Icons

import { library } from '@fortawesome/fontawesome-svg-core';
import { faUserSecret } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/vue-fontawesome';

library.add(faUserSecret);
Vue.component('font-awesome-icon', FontAwesomeIcon);

Vue.config.productionTip = false;

import SocketIO from  'socket.io-client'
import VueSocketIO from 'vue-socket.io'

Vue.use(new VueSocketIO({
    debug: true,
    connection: SocketIO('http://localhost:5000/')
  })
);


new Vue({
    el: '#app',
    router,
    render: h => h(App)
});
