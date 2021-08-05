const path = require('path');
const webpack = require('webpack');

module.exports = {
  entry: {
    'math-utils':'./src/index.ts',
  },
  devtool: 'source-map',
  node: { module: "empty", net: "empty", fs: "empty" },
  optimization:{
      minimize:false
  },
  output: {
    path: path.resolve(__dirname, 'dist'),
    filename: '[name].js',
    library: 'math-utils',
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