const path=require('path');
const webpack = require('webpack');

module.exports = {
    entry:'./src/index.ts',
    devtool: "inline-source-map",
    output:
    {
        filename:'bundle.js',
        path:path.resolve(__dirname,'dist'),
        /*hotUpdateChunkFilename: 'hot/hot-update.js',
        hotUpdateMainFilename: 'hot/hot-update.json',*/
        library: 'math',
        libraryTarget: 'var'
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