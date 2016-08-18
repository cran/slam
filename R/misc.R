
##
slam_options <-
local({
    options <- list(max_dense = 2^24)
    function(option, value) {
        if (missing(option)) return(options)
        if (missing(value))
            options[[option]]
        else
            options[[option]] <<- value
    }
})

