<template>

  <div id="app">
    
  <div class="app-body">
      <router-view/>
  </div>

  </div>
</template>

<script>

export default {
  name: 'App',
  data() {
    return {
      server: 'localhost:5000',
      navbar: true,
      active: false,
      active_info: false,
    };
  },
}
</script>

<style>
#app {
  background-color: #171E29;
}
#app {
  padding: 5rem;
  margin: 0;
  font-family: Consolas, monaco, monospace;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  width: 100%;
  min-height: 1000px; 
}
</style>
