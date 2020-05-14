module.exports = {
  runtimeCompiler: true,
  chainWebpack: config => {
        config.module
          .rule('vue')
          .use('vue-loader')
            .loader('vue-loader')
            .tap(options => {
              options.transformAssetUrls = {
                'vs-avatar': 'src'
              };
              return options
            })
    }
};