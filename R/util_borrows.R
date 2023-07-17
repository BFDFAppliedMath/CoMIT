#
# #FROM R.UTILS (https://rdrr.io/cran/R.utils/src/R/withTimeout.R)
# withTimeout <- function(expr, substitute=TRUE, envir=parent.frame(), timeout, cpu=timeout, elapsed=timeout, onTimeout=c("error", "warning", "silent"), ...) {
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   # Validate arguments
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   # Argument 'expr':
#   if (substitute) expr <- substitute(expr)
#
#   # Argument 'envir':
#   if (!is.environment(envir))
#     throw("Argument 'envir' is not a list: ", class(envir)[1L])
#
#   # Argument 'cpu' and 'elapsed':
#   cpu <- Arguments$getNumeric(cpu, range=c(0,Inf))
#   elapsed <- Arguments$getNumeric(elapsed, range=c(0,Inf))
#
#   # Argument 'onTimeout':
#   onTimeout <- match.arg(onTimeout)
#
#
#   setTimeLimit(cpu=cpu, elapsed=elapsed, transient=TRUE)
#   on.exit({
#     setTimeLimit(cpu=Inf, elapsed=Inf, transient=FALSE)
#   })
#
#   tryCatch({
#     eval(expr, envir=envir)
#   }, error = function(ex) {
#     msg <- ex$message
#     # Was it a timeout?
#     pattern <- gettext("reached elapsed time limit", "reached CPU time limit", domain="R")
#     pattern <- paste(pattern, collapse = "|")
#     if (regexpr(pattern, msg) != -1L) {
#       ex <- TimeoutException(msg, cpu=cpu, elapsed=elapsed)
#       if (onTimeout == "error") {
#         throw(ex)
#       } else if (onTimeout == "warning") {
#         warning(getMessage(ex))
#       } else if (onTimeout == "silent") {
#       }
#     } else {
#       # Rethrow error
#       throw(ex)
#     }
#   })
# }
