<template>
  <div>
  <div class="row">
    <div class="column-4">
    <gb-heading :tag="'h3'">NanoPath Live - Dashboard v0.1.0</gb-heading>
    <gb-divider id="summary_divider" color="purple" size="large"></gb-divider>
    <div id="summary">
        <gb-badge class="tag" color="turquoise">Online</gb-badge>
        <gb-badge class="tag" color="purple">DNA</gb-badge>
        <gb-badge class="tag" color="turquoise">Metagenomics</gb-badge>
        <p id="run_id">{{ run_id }}</p>
    </div>
    </div>
  </div>
  </div>
</template>

<style>
#summary {
    color: #A9C7DF;
}
#summary_divider {
    margin-top: 1rem;
    margin-bottom: 1rem;
}

#run_id {
    margin-top: 2rem;
}

.tag {
    margin-right: 1rem;
}


</style>


<script>

export default {
  name: 'Dashboard',
  data() {
    return {
        run_id: '7c254009-422b-46a2-8c69-d8a00d97363b'
    }
  },
  computed: {

  },
  sockets: {
    
  },
  methods: {

  },
  created() {

  },
};
</script>
