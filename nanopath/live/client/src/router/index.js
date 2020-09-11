import Vue from 'vue'
import Router from 'vue-router'
import Pathogen from '@/components/Pathogen'
import Signal from '@/components/Signal'


Vue.use(Router);

const router = new Router({
  routes: [
    {
        path: '*',
        redirect: '/public/404'
    },
    {
        path: '/',
        redirect: '/pathogen',
    },
    {
        path: '/pathogen',
        name: 'Pathogen',
        component: Pathogen,
    },
    {
      path: '/signal',
      name: 'Signal',
      component: Signal,
  },
  ],
  mode: 'history'
});

export default router;