const path = require('path');
const webpack = require('webpack');

module.exports = {
  entry: {
    'math-utils.min': './src/index.ts'
  },
  devtool:"source-map",
  node: { module: "empty", net: "empty", fs: "empty" },
  optimization:{
      minimize:true
  },
  output: {
    path: path.resolve(__dirname, 'dist'),
    filename: '[name].js',
    library: 'Exp',
    libraryTarget: 'var'
  },
  resolve:{
      extensions: ['.ts','.tsx','.js'],
  },
  module: {
      rules: [
          {
              test: /\.ts(x?)$/,
              exclude: /node_modules/,
              use: [
                  {
                      loader: "ts-loader"
                  }
              ]
          } 
      ]
  }
};