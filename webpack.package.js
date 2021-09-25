const path=require('path');
const webpack = require('webpack');

module.exports = {
    target:'node',
    entry:'./src/index.ts',
    devtool: "source-map",
    output: {
        filename:'math-utils-js.js',
        path:path.resolve(__dirname, './dist/package'),
        /*hotUpdateChunkFilename: 'hot/hot-update.js',
        hotUpdateMainFilename: 'hot/hot-update.json',*/
        library: 'MathUtilsJS',
        libraryTarget: 'commonjs'
    },
    plugins:
    [
	    //new webpack.HotModuleReplacementPlugin(),
    ],
    resolve: {
        extensions: [".ts", ".tsx",'.js']
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
    },
};