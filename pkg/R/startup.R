.onAttach <- function(lib, pkg) {
    ver <- read.dcf(file.path(lib,pkg,"DESCRIPTION"), "Version")
    msg <- sprintf("Development features for T-LoCoH\nVersion %s\nURL: http://tlocoh.r-forge.r-project.org/\nPlease send bug reports and feedback to: tlocoh@gmail.com", as.character(ver))
    packageStartupMessage(msg)
}
