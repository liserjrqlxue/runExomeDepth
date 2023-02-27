package main

import "github.com/liserjrqlxue/goUtil/stringsUtil"

type Mappability struct {
	chromosome      string
	start           int
	end             int
	meanMappability float64
}

func map2Map(item map[string]string) *Mappability {
	var m = &Mappability{
		chromosome:      item["chr"],
		start:           stringsUtil.Atoi(item["start"]),
		end:             stringsUtil.Atoi(item["end"]),
		meanMappability: stringsUtil.Atof(item["meanMappability"]),
	}
	return m
}
