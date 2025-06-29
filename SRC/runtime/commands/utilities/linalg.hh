
const char *linalg[] = {R"END(
namespace eval OpenSees {
namespace export verify

proc verify {cmd {value ""} {reference ""} {tolerance 1e-12} {about ""}} {
    if {$cmd == "error"} {
        set check [expr abs(($value - $reference)/$reference)]
        if {$check > $tolerance} {
        puts  "   \033\[31mFAIL\033\[0m: | $value - $reference | = $check > $tolerance"
        error "$about"
        } else {
          puts  "   \033\[32mPASS\033\[0m  "; # " $value   $reference $about"
        }

    } elseif {$cmd == "value"} {
        if {abs($value - $reference) > $tolerance} {
        set check [expr abs($value - $reference)]
        puts  "   \033\[31mFAIL\033\[0m($about): | $value - $reference | = $check > $tolerance"
        error "$about"
        } else {
        puts  "    \033\[32mPASS\033\[0m  "; # "$value   $reference $about"
        }
    } else {
     # value or  "about"
    puts "  $value"
    }
}
}

namespace import OpenSees::verify

proc range args {
  # https://wiki.tcl-lang.org/page/range
  foreach {start stop step} [switch -exact -- [llength $args] {
      1 {concat 0 $args 1}
      2 {concat   $args 1}
      3 {concat   $args  }
      default {error {wrong # of args: should be "range ?start? stop ?step?"}}
  }] break
  if {$step == 0} {error "cannot create a range when step == 0"}
  set range [list]
  while {$step > 0 ? $start < $stop : $stop < $start} {
      lappend range $start
      incr start $step
  }
  return $range
}

proc linspace {min  max num} {
    # https://opensource.apple.com/source/tcl/tcl-118.50.1/tcl_ext/tklib/tklib/examples/plotchart/rosenbrock.tcl.auto.htm

    #set opts {
    #  {min.arg 0.0 "the minimum value in the range"}
    #  {max.arg 1.0 "the maximum value in the range"}
    #  {step.arg 0.1 "the step to use"}
    #  {nb.arg 10 "the number of points in the range"}
    #  {method.arg "step" "the method to compute the values"}
    #}
    #set usage ": plotcontours \[options] ...\noptions:"
    #array set params [::cmdline::getoptions args $opts $usage]

    #set min $params(min)
    #set max $params(max)
    if {$min > $max} then {
      error "The minimum value is $min but the maximum value is $max"
    }
    set step [expr double($max-$min)/($num-1)]
    #switch -- $params(method) {
    #  "num" {
    #    set num $params(num)
    #    set step [expr {($max-$min)/($num-1)}]
    #  }
    #  "step" {
    #    set step $params(step)
    #  }
    #  default {
    #    error "Unknown method $params(method)"
    #  }
    #}
    set result {}
    set current $min
    while {1} {
      lappend result $current
      set current [expr {$current + $step}]
      if {$current > $max} then {
        break
      }
    }
    return $result
}
)END" ,
R"END(
# generator.tcl --
#
#       Iterators and generators via coroutines.
#
# Copyright (c) 2009 by Neil Madden <nem@cs.nott.ac.uk>
#
# See the file "license.terms" for information on Usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Tcl         8.6
package provide generator   0.2

namespace eval generator {
    namespace export {[a-z]*}
    namespace ensemble create

    # next generator varName ?varName ..? --
    #
    #   Fetch the next values from a generator, assigning them to variables. If
    #   the generator is exhausted any remaining variables are assigned the
    #   empty string.
    #
    proc next {generator args} {
        set items [takeList [llength $args] $generator]
        uplevel 1 [list lassign $items {*}$args]
    }

    proc takeList {n generator} {
        if {![exists $generator]} { return [list] }
        set ret [list]
        for {set i 0} {$i < $n} {incr i} {
            set item [$generator]
            if {[llength $item] == 0} { break }
            lappend ret [lindex $item 0]
        }
        return $ret
    }


    # foreach varSpec generator ?varSpec generator ...? body --
    #
    #       Iterate over the elements of one or more generator functions. Each
    #       generator is a command that yields successive elements. The syntax
    #       of this construct closely matches that of the built-in foreach
    #       command.
    #
    proc foreach args {
        if {[llength $args] < 3 || ([llength $args] % 2) != 1} {
            Usage "foreach varSpec generator ?varSpec generator ..? body"
        }
        set body [lindex $args end]
        set genSpec [lrange $args 0 end-1]
        set items [list]
        ::foreach {varList generator} $genSpec {
            lappend items [takeList [llength $varList] $generator]
        }

        try {
            # Keep going until all of the generators are exhausted
            # This is the foreach behaviour: empty strings are substituted for
            # exhausted generators.
            while 1 {
                set count 0
                ::foreach {varList generator} $genSpec item $items {
                    incr count [llength $item]
                    uplevel 1 [list lassign $item {*}$varList]
                }
                if {$count == 0} {
                    # All exhausted
                    break
                }
                try {
                    uplevel 1 $body
                } on continue {} {
                    # Continue processing
                } on return {result options} {
                    # increment -level to remove implementation details
                    dict incr options -level
                    return -options $options $result
                }
                set items [list]
                ::foreach {varList generator} $genSpec {
                    lappend items [takeList [llength $varList] $generator]
                }
            }
        } finally {
            # Ensure generators are all cleaned up
            ::foreach {_ generator} $genSpec {
                destroy $generator
            }
        }
        return
    }

    # finally cmd args.. --
    #
    #       Arranges for cmd to be called when the generator is destroyed. This
    #       can be used to perform cleanup in the event that a generator is
    #       terminated early.
    #
    proc finally args {
        set ns [uplevel 1 { namespace current }]
        trace add command [info coroutine] delete [list ::apply [list args $args $ns]]
    }
    proc exists generator { expr {[llength [info commands $generator]] != 0} }
    proc destroy args {
        ::foreach generator $args {
            if {[exists $generator]} { rename $generator "" }
        }
    }
    proc yield args {
        # Each argument is yielded individually as a separate value.
        ::foreach arg $args {
            ::yield [list $arg]
        }
    }

    proc define {name params body} {
        set name [Resolve 1 $name]
        set ns [namespace qualifiers $name]
        set lambda [list $params $body $ns]
        interp alias {} $name {} ::generator::spawn $lambda
        return $name
    }
                ##### PRIVATE METHODS #####

    proc spawn {lambda args} {
        set it [Gensym]
        coroutine $it ::generator::generate $it $lambda $args
    }

    proc generate {name lambda argList} {
        ::yield $name
        apply $lambda {*}$argList
    }

    proc Resolve {level name} {
        if {[string match ::* $name]} { return $name }
        if {[string is integer -strict $level] && $level >= 0} { incr level }
        set ns [uplevel $level { namespace current }]
        if {$ns eq "::"} { return ::$name }
        return $ns\::$name
    }

    proc All {p xs} {
        ::foreach x $xs { if {![{*}$p $x]} { return 0 } }
        return 1
    }

    proc Empty? xs { expr {[llength $xs] == 0} }


    variable Gensymid 0
    proc Gensym {} {
        variable Gensymid
        set prefix [namespace current]::generator
        while {1} {
            set name $prefix[incr Gensymid]
            if {[llength [info commands $name]] == 0} { break }
        }
        return $name
    }

    proc Usage msg {
        return -code error -level 2 "wrong # args: should be \"$msg\""
    }

    proc names {} {
        set pat {[0-9]*}
        return [info commands [namespace current]::generator$pat]
    }

        ##### STANDARD GENERATORS #####

    define map {f xs} {
        # Ensure underlying generator is cleaned up too
        finally destroy $xs
        foreach x $xs { yield [{*}$f $x] }
    }
    define filter {p xs} {
        finally destroy $xs
        foreach x $xs {
            if {[{*}$p $x]} { yield $x }
        }
    }
    proc reduce {f z xs} {
        foreach x $xs { set z [{*}$f $z $x] }
        return $z
    }
    proc foldl {f z xs} { reduce $f $z $xs }
    proc foldr {f z xs} {
        set ys [generator to list $xs]
        for {set i 0} {$i < [llength $ys]} {incr i} {
            set z [{*}$f [lindex $ys end-$i] $z]
        }
        return $z
    }
    define zipWith {f xs ys} {
        finally destroy $xs $ys
        foreach x $xs y $ys { yield [{*}$f $x $y] }
    }
    proc zip {xs ys} { zipWith list $xs $ys }

    proc all {p xs} {
        and [map $p $xs]
    }
    proc and xs {
        # foldl && true $xs
        # more efficient implementation (bail-out on first non-true element):
        foreach x $xs { if {!$x} { return 0 } }
        return 1
    }
    proc any {p xs} { reduce or 0 [map $p $xs] }
    define concat args {
        ::foreach xs $args { finally destroy $xs }
        ::foreach xs $args {
            foreach x $xs { yield $x }
        }
    }
    define concatMap {f xs} {
        concat {*}[map $f $xs]
    }
    proc drop {n xs} {
        takeList $n $xs
        return $xs
    }
    define dropWhile {p xs} {
        finally destroy $xs
        foreach x $xs {
            if {![{*}$p $x]} { yield $x; break }
        }
        foreach x $xs { yield $x }
    }
    proc contains {elem xs} {
        foreach x $xs { if {$x eq $elem} { return 1 } }
        return 0
    }

    proc foldl1 {f xs} { foldl $f [take 1 $xs] $xs }
    proc foldli {f z xs} {
        foreach x $xs { set z [{*}$f [incr i] $z $x] }
        return $z
    }
    proc foldri {f z xs} {
        set ys [to list $xs]
        for {set i [llength $ys]} {$i > 0} {incr i -1} {
            set z [{*}$f [incr j] [lindex $ys $i-1] $z]
        }
        return $z
    }
    proc head xs { take 1 $xs }
    proc tail xs { drop 1 $xs }
    proc last xs {
        foreach x $xs { }
        return $x
    }
    define init xs {
        finally destroy $xs
        set last [head $xs]
        foreach x $xs {
            yield $last
            set last $x
        }
    }
    define take {n xs} {
        finally destroy $xs
        foreach x $xs {
            if {[incr i] >= $n} { break }
            yield $x
        }
    }
    define iterate {f x} {
        while 1 {
            yield $x
            set x [{*}$f $x]
        }
    }

    proc Count {x y} { incr x }
    proc length xs { foldl Count 0 $xs }

    proc or {p xs} {
        foreach x $xs { if {[{*}$p $x]} { return 1 } }
        return 0
    }

    proc product xs { foldl ::tcl::mathop::* 1 $xs }
    define repeat {n args} {
        for {set i 0} {$i < $n} {incr i} {
            yield {*}$args
        }
    }
    proc sum xs { foldl ::tcl::mathop::+ 0 $xs }

    define takeWhile {p xs} {
        finally destroy $xs
        foreach x $xs {
            if {[{*}$p $x]} {
                yield $x
            } else {
                break
            }
        }
    }

    define splitWhen {p xs} {
        finally destroy $xs
        set token [list]
        foreach x $xs {
            if {[{*}$p $x]} {
                yield $token
                set token [list]
            } else {
                lappend token $x
            }
        }
        if {[llength $token]} { yield $token }
    }

    define scanl {f z xs} {
        finally destroy $xs
        yield $z
        foreach x $xs {
            set z [{*}$f $z $x]
            yield $z
        }
    }

    # from ?list|dict? xs
    # Converts a list or dictionary into a generator: over elements or key/value
    # pairs.
    namespace eval from {
        namespace export list dict string
        namespace ensemble create

        generator define list xs {
            foreach x $xs { generator yield $x }
        }

        generator define dict d {
            ::dict for {k v} $d { generator yield $k $v }
        }

        generator define string s {
            foreach c [split $s ""] { generator yield $c }
        }
    }

    # to ?list|dict? g
    # Converts a generator into a list or dictionary by extracting all elements.
    # Dictionaries are created by assuming the generator returns a pair of
    # values per element, and using these as the key and value.
    namespace eval to {
        namespace export list dict string
        namespace ensemble create

        proc list g {
            set xs [::list]
            generator foreach x $g { lappend xs $x }
            return $xs
        }

        proc dict g {
            set d [dict create]
            generator foreach {k v} $g { dict set d $k $v }
            return $d
        }

        proc string g {
            set s ""
            generator foreach c $g { append s $c }
            return $s
        }
    }
    # The conversion functions should follow these identity laws:
    # [to list [from list $xs]] == $xs
    # [to dict [from dict $xs]] == $xs
}
)END"
};

