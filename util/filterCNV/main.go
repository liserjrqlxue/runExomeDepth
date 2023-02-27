package main

import (
	"bufio"
	"embed"
	"flag"
	"github.com/liserjrqlxue/goUtil/fmtUtil"
	"github.com/liserjrqlxue/goUtil/math"
	"github.com/liserjrqlxue/goUtil/osUtil"
	"github.com/liserjrqlxue/goUtil/scannerUtil"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/liserjrqlxue/goUtil/stringsUtil"
	"github.com/liserjrqlxue/goUtil/textUtil"
	"log"
)

var (
	in = flag.String(
		"in",
		"",
		"input cnv",
	)
	out = flag.String(
		"out",
		"",
		"output cnv",
	)
	gender = flag.String(
		"gender",
		"",
		"gender",
	)
)

//go:embed map.txt
var emFS embed.FS

func main() {
	flag.Parse()
	if *in == "" || *out == "" {
		log.Fatal("-in/-out required!")
	}

	var meanMappability = make(map[string][]*Mappability)
	var mapTxt, err = emFS.Open("map.txt")
	simpleUtil.CheckErr(err)
	defer simpleUtil.DeferClose(mapTxt)
	var scan = bufio.NewScanner(mapTxt)
	var maps, _ = scannerUtil.Scanner2MapArray(scan, "\t", nil)

	for _, item := range maps {
		var m = map2Map(item)
		var ms = meanMappability[m.chromosome]
		ms = append(ms, m)
		meanMappability[m.chromosome] = ms
	}

	var cnv, title = textUtil.File2MapArray(*in, "\t", nil)

	var output = osUtil.Create(*out)
	defer simpleUtil.DeferClose(output)

	fmtUtil.FprintStringArray(output, title, "\t")

	switch *gender {
	case "M":
		for _, m := range cnv {
			var (
				ratio     = stringsUtil.Atof(m["reads.ratio"])
				bf        = stringsUtil.Atof(m["BF"])
				chr       = m["chromosome"]
				start     = stringsUtil.Atoi(m["start"])
				end       = stringsUtil.Atoi(m["end"])
				meanMs    []float64
				meanScore = 1.0
			)
			if ratio > 0.55 && ratio < 1.35 {
				continue
			}
			if bf < 10 && ratio > 0.5 && ratio < 1.5 {
				continue
			}
			for _, meanM := range meanMappability[chr] {
				if meanM.start <= end && meanM.end >= start {
					meanMs = append(meanMs, meanM.meanMappability)
				}
			}
			if len(meanMs) > 0 {
				meanScore = math.Mean(meanMs)
			}
			if bf < 55 && meanScore < 0.8 {
				continue
			}

			var line []string
			for _, s := range title {
				line = append(line, m[s])
			}
			fmtUtil.FprintStringArray(output, line, "\t")
		}
	default:
		for _, m := range cnv {
			var line []string
			for _, s := range title {
				line = append(line, m[s])
			}
			fmtUtil.FprintStringArray(output, line, "\t")
		}
	}
}
