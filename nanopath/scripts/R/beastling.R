
process_log <- function(
    log_file=NULL, burnin=0.1, hpd=0.95, most_recent=2019, contemporary=FALSE, sample_slice=FALSE
) {

    if (contemporary){
        # No null intervals in Contemporary Rho sampling prop (BEAST not happy)
        sampling_prop_label = "rho_BDSKY_Contemp"  # unless it is Birth-Death Skyline Contemporary
        become_uninfectious_rate_label = "becomeUninfectiousRate_BDSKY_Contemp"
        origin_label = "origin_BDSKY_Contemp"
        reproduction_number_label = "reproductiveNumber_BDSKY_Contemp"
    } else {
        # Otherwise in Serial sampling proportion get all parameters

        origin_label = "origin_BDSKY_Serial"
        reproduction_number_label = "reproductiveNumber_BDSKY_Serial"
        become_uninfectious_rate_label = "becomeUninfectiousRate_BDSKY_Serial"

        if (sample_slice) {
            sampling_prop_label = "samplingProportion_BDSKY_Serial.2" # post-sample slice, pre-sample zero fixed
        } else {
            sampling_prop_label = "samplingProportion_BDSKY_Serial"
        }

    }

    run <- readLogfile(log_file, burnin=burnin)

    r0 <- getSkylineSubset(run,  "reproductiveNumber")
    r0_hpd <- getMatrixHPD(r0, alpha=1-hpd)

    # This is conditioning the time scale on the median tree height
    median_age <- median(run$TreeHeight)
    median_origin <- most_recent-median_age # e.g: 1986.82

    # Skyline grid: length of years until median posterior MRCA
    timegrid <- 1:ceiling(median_age)  # round up to capture all time points
    age_samples <- run$TreeHeight

    r0_gridded <- gridSkyline(r0, age_samples, timegrid)
    r0_gridded_hpd <- getMatrixHPD(r0_gridded)

    # R0 Skyline data (can be un-sliced)
    hpd_data<- t(r0_gridded_hpd)
    hpd_data <- as.data.frame(
        cbind(hpd_data, as.integer(median_origin)+timegrid)
    )
    # Intervals
    colnames(hpd_data) <- c('Lower', 'Mean', 'Upper', 'Year')

    print(
        paste(
            become_uninfectious_rate_label,
            origin_label, sampling_prop_label, most_recent
        )
    )

    # Posterior data
    df = data.frame(
        TreeHeight=most_recent-run$TreeHeight,
        ClockRate=run$clockRate,
        InfectiousPeriod=1/run[[become_uninfectious_rate_label]],
        Origin=most_recent-run[[origin_label]],
        SamplingProportion=run[[sampling_prop_label]]
    )

    if (ncol(r0) == 1){
        df$ReproductionNumber <- run[[reproduction_number_label]]
    }

    results = list(
        "posterior_data"=df,
        "reproduction_number_skyline"=hpd_data
    )

    return(results)
}

get_hpd_title <- function(data, prefix='Origin', dig=3, science=FALSE, suffix=""){

    hpd <- getHPD(data)

    if (science){
        lower <- scales::scientific(hpd[1])
        median <-  scales::scientific(hpd[2])
        upper <- scales::scientific(hpd[3])
    } else {
        lower <- round(hpd[1], dig)
        median <-  round(hpd[2], dig)
        upper <- round(hpd[3], dig)
    }

    if (var(data) == 0){
        hpd_title <- "fixed"
    } else {
        hpd_title <- paste0("95% HPD: ", lower, " - ", upper)
    }
    title <- paste0(prefix, "\n", median, " ", suffix, " (", hpd_title, ")\n")

    return(title)
}

plot_single_log <- function(posterior_data, reproduction_number_skyline, oldest, plot_path, plot_name, name=NULL, hpds=NULL){


    p1_title <- get_hpd_title(posterior_data$Origin, prefix='Origin', dig=2)
    p1 <- ggplot(posterior_data, aes(Origin)) + ggtitle(p1_title) +
        geom_density(alpha=1, fill="#2171b5") + xlab("") + ylab("") +
        theme_cowplot() + theme(axis.title.y = element_blank(),
                                axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                axis.line.y = element_blank(), panel.grid.minor = element_blank(),
                                plot.title = element_text(size=12, face="plain")
        )

    p2_title <- get_hpd_title(posterior_data$TreeHeight, prefix='MRCA', dig=2)
    p2 <- ggplot(posterior_data, aes(TreeHeight)) + ggtitle(p2_title) +
        geom_density(alpha=1, fill="#54278f") + xlab("") + ylab("") +
        theme_cowplot() + theme(axis.title.y = element_blank(),
                                axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                axis.line.y = element_blank(), panel.grid.minor = element_blank(),
                                plot.title = element_text(size=12, face="plain")
        )


    p3_title <- get_hpd_title(posterior_data$ClockRate, prefix='Clock rate', science=TRUE)
    p3 <- ggplot(posterior_data, aes(ClockRate)) + ggtitle(p3_title) + scale_x_continuous(labels = scales::scientific) +
        geom_density(alpha=1, fill="#1d91c0") + xlab("") + ylab("") +
        theme_cowplot() + theme(axis.title.y = element_blank(),
                                axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                axis.line.y = element_blank(), panel.grid.minor = element_blank() ,
                                plot.title = element_text(size=12, face="plain")
        )


    xupper <- max(posterior_data$InfectiousPeriod)
    if ( xupper <= 2) {
        xstep <- 0.2
    } else if (xupper < 5){
        xstep <- 0.5
    } else if (xupper < 10){
        xstep <- 10
    } else if (xupper > 10 && xupper <= 50) {
        xstep <- 5
    } else {
        xstep <- 10
    }


    p4_title <- get_hpd_title(posterior_data$InfectiousPeriod, prefix='Infectious period', suffix='years')
    p4 <- ggplot(posterior_data, aes(InfectiousPeriod)) + ggtitle(p4_title) +
        scale_x_continuous(breaks=seq(0, ceiling(max(posterior_data$InfectiousPeriod)), xstep)) +
        geom_density(alpha=1, fill="#238b45") + xlab("") + ylab("years") +
        theme_cowplot() + theme(axis.title.y = element_blank(),
                                axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                axis.line.y = element_blank(), panel.grid.minor = element_blank(),
                                plot.title = element_text(size=12, face="plain")
        )

    xupper <- max(posterior_data$SamplingProportion)
    if ( xupper <= 0.1){
        xstep <- 0.01
        xmax <- 0.1
    } else {
        xstep <- 0.1
        xmax <- 1
    }

    p5_title <- get_hpd_title(posterior_data$SamplingProportion, prefix='Sampling proportion')
    p5 <- ggplot(posterior_data, aes(SamplingProportion)) + ggtitle(p5_title) +
        scale_x_continuous(breaks=seq(0, xmax, xstep)) +
        geom_density(alpha=1, fill="#fd8d3c") + xlab("") + ylab("") + theme_cowplot() +
        theme_cowplot() + theme(axis.title.y = element_blank(),
                                axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                axis.line.y = element_blank(), panel.grid.minor = element_blank(),
                                plot.title = element_text(size=12, face="plain")
        )

    # Reproduction number interval plot
    if (is.null(reproduction_number_skyline)) {
        # Standard posterior plot:
        r0_title <- get_hpd_title(posterior_data$ReproductionNumber, prefix='Reproduction number')
        r0_plot <- ggplot(posterior_data, aes(ReproductionNumber)) + ggtitle(r0_title) +
            scale_x_continuous(breaks=seq(0, ceiling(max(posterior_data$ReproductionNumber)), 1)) +
            geom_density(alpha=1, fill="#f03b20") + xlab("") + ylab("") + theme_cowplot() +
            theme(
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.y = element_blank(),
                legend.title = element_blank(),
                panel.grid.minor = element_blank(),
                plot.title = element_text(size=12, face="plain")
            )
    } else {
        # Interval plot over time:

        # r0_title <- get_hpd_title(posterior_data$ReproductionNumber, prefix='Reproduction number')
        # print(r0_title)

        ymax <- ceiling(max(reproduction_number_skyline$Upper))
        r0_plot = ggplot(reproduction_number_skyline, aes(x=Year, y=Mean)) + ylab("R0") +
            geom_ribbon(aes(ymin=Lower, ymax=Upper), fill="#d3d3d3") + ylab("") + xlab("") + ylim(0, ymax) +
            geom_line(color='black') + theme_classic() + ggtitle("R0\n") + geom_hline(yintercept=1, linetype="dotted", size=0.5)

        if (!is.null(oldest)){  # oldest just adds the line here
            r0_plot = r0_plot + geom_vline(
                xintercept = oldest, linetype="dotted", size=0.5
            )
        }

        if (!(is.null(hpds))){
            hpd_data <- hpds[hpds$name == name,]

            xmin <- as.numeric(hpd_data$xmin)
            xmax <- as.numeric(hpd_data$xmax)
            color <- as.character(hpd_data$color)

            for (i in seq(nrow(hpd_data))){
                r0_plot <- r0_plot + geom_rect(xmin=xmin[i], xmax=xmax[i], ymin=0, ymax=ymax, alpha=0.01, fill=color[i])

            }
        }
    }

    po <- cowplot::plot_grid(p2, p1, p3, p4, p5, r0_plot)
    cowplot::save_plot(paste0(plot_path, "/", plot_name, ".pdf"), po, base_asp=2.1, base_height=7.7)
    cowplot::save_plot(paste0(plot_path, "/", plot_name,".svg"), po, base_asp=2.1, base_height=7.7)

}


# Summary functions

plot_posterior_single <- function(log_file, oldest, most_recent, burnin=0.1, hpd=0.95, sample_slice=TRUE){

    # Summary function to plot the posterior distributions of a single log file
    # of a BDSKY Serial or Contemporary model run in BEAST2

    if (oldest == most_recent){
        contemporary = TRUE
    } else {
        contemporary = FALSE
    }

    results <- process_log(
        log_file=log_file, burnin=burnin, hpd=hpd, sample_slice=sampel_slice,
        most_recent=most_recent, contemporary=contemporary
    )

    plot_single_log(
        posterior_data = results$posterior_data,
        reproduction_number_skyline = NULL,  # posterior plot, no skyline plot
        oldest=oldest_bdss, name=tools::file_path_sans_ext(log_file)
    )

}

read_log_files <- function(
    data, path=NULL, file=NULL, burnin=01., hpd=0.95,
    extract_label_function=NULL, extract_dimension_function=NULL,
    logpattern=".log", sample_slice=TRUE
){

    if (!is.null(path)){
        files = list.files(path=path, pattern=logpattern, full.names=TRUE)
        files = stringr::str_sort(files, numeric=TRUE)
    } else {
        files = c(file)
    }

    processed <- lapply(files, function(f){

        # Use a data frame to specify data per label

        file_name <- basename(f)
        fname <- tools::file_path_sans_ext(file_name)

        if (!is.null(extract_label_function)){
            name <- extract_label_function(fname)
        }
        #print("Name")
        print(name)

        if (!is.null(extract_dimension_function)){
            d <- extract_dimension_function(fname)
        } else {
            d <- NULL
        }

        #print(d)
        #print(data)
        #print(name)

        most_recent = data[name, "last_sample_year"]
        oldest = data[name, "first_sample_year"]

        #print(most_recent)
        #print(oldest)

        if (most_recent == oldest){
            contemporary = TRUE
        } else {
            contemporary = FALSE
        }

        results <- process_log(
            log_file=f, burnin=burnin, hpd=hpd, sample_slice=sample_slice,
            most_recent=most_recent, contemporary=contemporary
        )


        if (!is.null(d)){
            name <- paste0(name, "_", d)
        }

        results$posterior_data$run_name <- rep(
            name, length(results$posterior_data$TreeHeight)
        )

        results$reproduction_number_skyline$run_name <- rep(
            name, length(results$reproduction_number_skyline$Year)
        )

        results$meta_data <- data.frame(
            name=name, oldest=oldest, most_recent=most_recent, file=f
        )


        return(results)

    })

    posteriors <- lapply(processed, function(x) { return(x$posterior_data) })
    df1 <- do.call(rbind, posteriors)

    df1$run_name <- factor(
        df1$run_name, levels=stringr::str_sort(
            levels(factor(df1$run_name)), numeric=TRUE
        )
    )


    re <- lapply(processed, function(x) { return(x$reproduction_number_skyline) })
    df2 <- do.call(rbind, re)

    df2$run_name <- factor(
        df2$run_name, levels=stringr::str_sort(
            levels(factor(df2$run_name)), numeric=TRUE
        )
    )

    metas <- lapply(processed, function(x) { return(x$meta_data) })
    df3 <- do.call(rbind, metas)

    results = list(
        "posterior"=df1,
        "re"=df2,
        "meta"=df3
    )

    return(results)


}

# Plotting functions

create_sensitivity_plots <- function(df, mrca_limits, origin_limits, plot_path, plot_name){

    p1 <- ggplot(df, aes(x=TreeHeight, fill=sensitivity)) + ggtitle("MRCA\n") +
        geom_density(alpha=.7) + scale_fill_manual(values=c("prior"="#D3D3D3", "posterior"="#54278f")) + xlab("") + ylab("") +
        theme_cowplot() + theme(
            panel.grid.minor = element_blank(),
            legend.position="none"
        )

    if (!is.null(mrca_limits)){
        p1 <- p1 + xlim(mrca_limits)
    } else {
        p1 <- p1 + xlim(floor(min(df$TreeHeight)), ceiling(max(df$TreeHeight))+1)
    }

    p2 <- ggplot(df, aes(x=Origin, fill=sensitivity)) + ggtitle("Origin\n") +
        geom_density(alpha=.7) + scale_fill_manual(values=c("prior"="#D3D3D3", "posterior"="#08519c")) + xlab("") + ylab("") +
        theme_cowplot() + theme(
            panel.grid.minor = element_blank(),
            legend.position="none"
        )

    if (!is.null(origin_limits)){
        p2 <- p2 + xlim(origin_limits)
    } else {
        p2 <- p2 + xlim(floor(min(df$Origin)), ceiling(max(df$Origin))+1)
    }

    if (var(df$ClockRate) == 0){
        clock_rate <- unique(df$ClockRate)[1]
        min_clock <- clock_rate - 0.5*clock_rate
        max_clock <- clock_rate + 0.5*clock_rate
        p3 <- ggplot() + ggtitle("Clock rate (fixed)\n") + scale_x_continuous(labels = scales::scientific) +
            xlab("") + ylab("") + geom_vline(xintercept = clock_rate, linetype="dotted", size=1) +
            xlim(min_clock, max_clock)  +
            theme_cowplot() + theme(
                panel.grid.minor = element_blank(),
                legend.position="none"
            )
    } else {
        p3 <- ggplot(df, aes(x=ClockRate, fill=sensitivity)) + ggtitle("Clock rate (SNP / site / year)\n") + scale_x_continuous(labels = scales::scientific) +
            geom_density(alpha=.7) + scale_fill_manual(values=c("prior"="#D3D3D3", "posterior"="#7bccc4")) + xlab("") + ylab("") +
            theme_cowplot() + theme(
                panel.grid.minor = element_blank(),
                legend.position="none"
            )
    }


    xupper <- max(df$InfectiousPeriod)
    if ( xupper <= 2) {
        xstep <- 0.2
    } else if (xupper < 5){
        xstep <- 0.5
    } else if (xupper < 10){
        xstep <- 10
    } else if (xupper > 10 && xupper <= 50) {
        xstep <- 5
    } else {
        xstep <- 10
    }


    p4 <- ggplot(df, aes(x=InfectiousPeriod, fill=sensitivity)) + ggtitle("Infectious period (years)\n") +
        geom_density(alpha=.7) + scale_fill_manual(values=c("prior"="#D3D3D3", "posterior"="#006d2c")) + xlab("") + ylab("") +
        scale_x_continuous(breaks=seq(0, ceiling(max(df$InfectiousPeriod)), xstep)) +
        theme_cowplot() + theme(
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(size=12),
            legend.position="none"
        )

    xupper <- max(df$SamplingProportion)
    if ( xupper <= 0.1){
        xstep <- 0.01
        xmax <- 0.1
    } else {
        xstep <- 0.1
        xmax <- 1
    }

    p5 <- ggplot(df, aes(x=SamplingProportion, fill=sensitivity)) + ggtitle("Sampling proportion\n") +
        geom_density(alpha=.7) + scale_fill_manual(values=c("prior"="#D3D3D3", "posterior"="#f16913")) + xlab("") + ylab("") +
        scale_x_continuous(breaks=seq(0, xupper, xstep)) +
        theme_cowplot() + theme(
            panel.grid.minor = element_blank(),
            legend.position="none"
        )

    if ("ReproductionNumber" %in% colnames(df)){

        p6 <- ggplot(df, aes(x=ReproductionNumber, fill=sensitivity)) + ggtitle("Reproduction number\n") +
            geom_density(alpha=.7) + scale_fill_manual(values=c("prior"="#D3D3D3", "posterior"="#e31a1c")) + xlab("") + ylab("") +
            scale_x_continuous(breaks=seq(0, ceiling(max(df$ReproductionNumber)), 1)) +
            geom_vline(xintercept = 1, linetype="dotted", size=0.5) +
            theme_cowplot() + theme(
                panel.grid.minor = element_blank(),
                legend.position="none"
            )
    } else {
        p6 <- NULL
    }

    po <- cowplot::plot_grid(p2, p1, p3, p4, p5, p6)
    cowplot::save_plot(paste0(plot_path, "/", plot_name, ".pdf"), po, base_asp=2.1, base_height=7.7)
    cowplot::save_plot(paste0(plot_path, "/", plot_name, ".png"), po, base_asp=2.1, base_height=7.7)


}

create_ridge_plots <- function(df, plot_path, plot_name, infectious_limit=NULL){

    p1 <- ggplot(df, aes(x=TreeHeight, y=run_name, fill=run_name)) + ggtitle("MRCA\n") +
        geom_density_ridges(scale=1.5, alpha=1) + scale_fill_brewer(palette = "Purples") + xlab("") + ylab("") +
        theme_cowplot() + theme(
            panel.grid.minor = element_blank(),
            legend.position="none"
        )

    p2 <- ggplot(df, aes(x=Origin, y=run_name, fill=run_name)) + ggtitle("Origin\n") +
        geom_density_ridges(scale=1.5, alpha=1) + scale_fill_brewer(palette = "Blues") + xlab("") + ylab("") +
        theme_cowplot() + theme(
            panel.grid.minor = element_blank(),
            legend.position="none"
        )

    if (var(df$ClockRate) == 0){
        clock_rate <- unique(df$ClockRate)[1]
        min_clock <- clock_rate - 0.5*clock_rate
        max_clock <- clock_rate + 0.5*clock_rate
        p3 <- ggplot() + ggtitle("Clock rate (fixed)\n") + scale_x_continuous(labels = scales::scientific) +
            xlab("") + ylab("") + geom_vline(xintercept = clock_rate, linetype="dotted", size=1) +
            xlim(min_clock, max_clock)  +
            theme_cowplot() + theme(
                panel.grid.minor = element_blank(),
                legend.position="none"
            )
    } else {
        p3 <- ggplot(df, aes(x=ClockRate, y=run_name, fill=run_name)) + ggtitle("Clock rate (SNP / site / year)\n") + scale_x_continuous(labels = scales::scientific) +
            geom_density_ridges(scale=1.5, alpha=1) + scale_fill_brewer(palette = "YlGnBu") + xlab("") + ylab("") +
            theme_cowplot() + theme(
                panel.grid.minor = element_blank(),
                legend.position="none"
            )
    }

    xupper <- ceiling(max(df$InfectiousPeriod))
    if ( xupper <= 2){
        xstep <- 0.2
    } else if (xupper < 5){
        xstep <- 0.5
    } else if (xupper < 10){
        xstep <- 10
    } else if (xupper > 10 && xupper <= 50) {
        xstep <- 5
    } else {
        xstep <- 10
    }
    if (!is.null(infectious_limit)){
        xmax <- infectious_limit
    } else {
        xmax <- ceiling(max(df$InfectiousPeriod))
    }
    p4 <- ggplot(df, aes(x=InfectiousPeriod, y=run_name, fill=run_name)) + ggtitle("Infectious period (years)\n") +
        geom_density_ridges(scale=1.5, alpha=1) + scale_fill_brewer(palette = "Greens") + xlab("") + ylab("") +
        scale_x_continuous(breaks=seq(0, xmax, xstep)) + xlim(0, xmax) +
        theme_cowplot() + theme(
            panel.grid.minor = element_blank(),
            legend.position="none"
        )



    xupper <- max(df$SamplingProportion)
    if ( xupper <= 0.1){
        xstep <- 0.01
        xmax <- 0.1
    } else {
        xstep <- 0.1
        xmax <- 1
    }

    p5 <- ggplot(df, aes(x=SamplingProportion, y=run_name, fill=run_name)) + ggtitle("Sampling proportion\n") +
        geom_density_ridges(scale=1.5, alpha=1) + scale_fill_brewer(palette = "Oranges") + xlab("") + ylab("") +
        scale_x_continuous(breaks=seq(0, xmax, xstep)) +
        theme_cowplot() + theme(
            panel.grid.minor = element_blank(),
            legend.position="none"
        )

    if ("ReproductionNumber" %in% colnames(df)){

        p6 <- ggplot(df, aes(x=ReproductionNumber, y=run_name, fill=run_name)) + ggtitle("Reproduction number\n") +
            geom_density_ridges(scale=1.5, alpha=1) + scale_fill_brewer(palette = "YlOrRd") + xlab("") + ylab("") +
            scale_x_continuous(breaks=seq(0, ceiling(max(df$ReproductionNumber)), 1)) +
            geom_vline(xintercept = 1, linetype="dotted", size=0.5) +
            theme_cowplot() + theme(
                panel.grid.minor = element_blank(),
                legend.position="none"
            )

    } else {
        p6 <- NULL
    }

    po <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6)

    cowplot::save_plot(paste0(plot_path, "/", plot_name, ".pdf"), po, base_asp=2.1, base_height=7.7)
    cowplot::save_plot(paste0(plot_path, "/", plot_name, ".svg"), po, base_asp=2.1, base_height=7.7)
}


plot_posterior_multiple <- function(df_posterior, plot_path, plot_name, infectious_limit=NULL){

    ridges <- create_ridge_plots(df=df_posterior, plot_path=plot_path, plot_name=plot_name, infectious_limit=infectious_limit)

}

plot_prior_sensitivity <- function(
    df_prior, df_posterior, plot_path,
    mrca_limits=NULL, origin_limits=NULL
){

    run_names <- levels(df_prior$run_name)
    # Sensitivity plot loop
    for (run in run_names){

        priors <- df_prior[df_prior$run_name == run, ]
        posteriors <- df_posterior[df_posterior$run_name == run, ]

        priors$sensitivity <- rep("prior", nrow(priors))
        posteriors$sensitivity <- rep("posterior", nrow(posteriors))

        df <- rbind(priors, posteriors)

        create_sensitivity_plots(
            df,
            mrca_limits=mrca_limits,
            origin_limits=origin_limits,
            plot_path=plot_path,
            plot_name=run
        )
    }
}




