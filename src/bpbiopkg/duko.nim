import times
import duktape/js

template len*(ctx: DTContext): int =
    ## return the size of the duktape stack
    ctx.duk_get_top()

template pop*(ctx:DTContext) =
    ## remove the top item from the duktape stack
    ctx.duk_pop()

proc `[]=`*(ctx: DTContext, key:string, value: SomeOrdinal) {.inline.} =
    ## set a global value in the context
    ctx.duk_push_int(value.duk_int_t)
    discard ctx.duk_put_global_string(key)

proc `[]=`*(ctx: DTContext, key: string, value: SomeFloat) {.inline.} =
    ## set a global value in the context
    ctx.duk_push_number(value.duk_double_t)
    discard ctx.duk_put_global_string(key)

proc `[]=`*(ctx: DTContext, key: string, values: seq[SomeNumber]) {.inline.} =
  var idx = ctx.duk_push_array()
  for i, v in values:
    ctx.duk_push_number(v.duk_double_t)
    discard ctx.duk_put_prop_index(idx, i.duk_uarridx_t)
  discard ctx.duk_put_global_string(key)

proc `[]`*(ctx: DTContext, key:string): float {.inline.} =
    if ctx.duk_get_global_string(key) == 0:
        raise newException(KeyError, "couldn't find key:" & key)
    result = ctx.duk_get_number(-1).float

proc check*(ctx: DTContext, expression: string): bool {.inline.} =
    ## evaluate the expression in the current context
    ctx.duk_eval_string(expression)
    result = ctx.duk_get_boolean(-1) == 1
    ctx.pop()

type Duko* = object
    ctx*: DTContext
    name*: string
    vptr: pointer

proc newObject*(ctx:DTContext, name: string): Duko =
  ## create a new object.
  result = Duko(ctx: ctx, name: name)
  discard ctx.duk_push_object()
  result.vptr = ctx.duk_get_heapptr(-1)
  doAssert result.ctx.duk_put_global_string(name) == 1

proc `[]=`*(o: Duko, key:string, value: SomeFloat) {.inline.} =
    ## set the property at key to a value
    var idx = o.ctx.duk_push_heapptr(o.vptr)
    o.ctx.duk_push_number(value.duk_double_t)
    doAssert o.ctx.duk_put_prop_string(idx, key) == 1
    o.ctx.pop()

proc clear*(o: Duko) {.inline.} =
    # TODO
    discard

proc `[]=`*(o: Duko, key:string, value: SomeInteger) {.inline.} =
    ## set the property at key to a value
    var idx = o.ctx.duk_push_heapptr(o.vptr)
    o.ctx.duk_push_int(value.duk_int_t)
    doAssert o.ctx.duk_put_prop_string(idx, key) == 1
    o.ctx.pop()

proc `[]=`*(o: Duko, key: string, values: seq[SomeNumber]) {.inline.} =
  var idx = o.ctx.duk_push_heapptr(o.vptr)
  var arr_idx = o.ctx.duk_push_array()
  for i, v in values:
    o.ctx.duk_push_number(v.duk_double_t)
    discard o.ctx.duk_put_prop_index(arr_idx, i.duk_uarridx_t)
  doAssert o.ctx.duk_put_prop_string(idx, key) == 1
  o.ctx.pop()

proc `[]`*(o: Duko, key:string): float {.inline.} =
    var idx = o.ctx.duk_push_heapptr(o.vptr)
    discard o.ctx.duk_push_string(key)
    doAssert o.ctx.duk_get_prop(idx) == 1
    result = o.ctx.duk_get_number(-1).float
    o.ctx.duk_pop_n(2)

when isMainModule:
  import unittest

  suite "duktaper":
    test "usage":
      var ctx = duk_create_heap_default()
      var kid = ctx.newObject("kid")
      kid["xyz"] = 22.31
      kid["asdf"] = 33.333
      check kid["xyz"] == 22.31
      check kid["asdf"] == 33.333
      kid["asdf"] = 33
      check kid["asdf"] == 33

      check ctx.len == 0
      ctx["somevar"] = 22.33
      check ctx.len == 0
      check ctx["somevar"] == 22.33
      ctx.duk_destroy_heap();


    test "array":
      var ctx = duk_create_heap_default()
      ctx["arr"] = @[22, 33]
      ctx.duk_eval_string("arr[0]")
      check ctx.duk_get_number(-1) == 22

      var dad = ctx.newObject("dad")
      dad["arr"] = @[55.2, 66.6, 22.3]
      ctx.duk_eval_string("dad.arr[1]")
      check ctx.duk_get_number(-1) == 66.6

    test "speed":
      var ctx = duk_create_heap_default()
      var kid = ctx.newObject("kid")

      var t = cpuTime()
      var tries = 100_000

      var success = 0
      for i in 0..tries:
          kid["dp"] = i
          ctx["mom"] = i.float
          ctx["dad"] = 23
          ctx["proband"] = i.float

          if ctx.check("kid.dp < 500 && dad > 21 && mom == kid.dp"):
            success.inc

      check success == 500
      echo (cpuTime() - t) / (tries / 1000000), " seconds per million evaluations"
      ctx.duk_destroy_heap();
