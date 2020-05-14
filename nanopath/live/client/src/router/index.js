import Vue from 'vue'
import Router from 'vue-router'
import Dashboard from '@/components/Dashboard'

Vue.use(Router);

const router = new Router({
  routes: [
    {
        path: '*',
        redirect: '/public/404'
    },
    {
        path: '/',
        redirect: '/dashboard',
    },
    {
        path: '/dashboard',
        name: 'Dashboard',
        component: Dashboard,
    },
  ],
  mode: 'history'
});

export default router;