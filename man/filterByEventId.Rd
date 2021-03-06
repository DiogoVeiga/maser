% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterEvents.R
\name{filterByEventId}
\alias{filterByEventId}
\title{Filter splicing events based on event identifier and type.}
\usage{
filterByEventId(events, event_id, type = c("A3SS", "A5SS", "SE", "RI", "MXE"))
}
\arguments{
\item{events}{a maser object.}

\item{event_id}{numeric vector of event identifiers.}

\item{type}{character indicating splice type. Possible values are
\code{c("A3SS", "A5SS", "SE", "RI", "MXE")}.}
}
\value{
a maser object.
}
\description{
Filter splicing events based on event identifier and type.
}
\examples{
path <- system.file("extdata", file.path("MATS_output"), package = "maser")
hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
filterByEventId(hypoxia, 33208, "SE")
}
