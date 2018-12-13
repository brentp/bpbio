import times
import duktape/js
export js

type Duko* = object
    ctx*: DTContext
    name*: string
    vptr: pointer

type Dukexpr* = object
    ## a compiled expression
    expr*: string
    ctx*: DTContext
    vptr: pointer

converter toBool*(d: duk_bool_t): bool {.inline.} =
    ## automatically converts the duk_bool_t to the appropriate true/false value in nim.
    return d == 1

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
  ## set a global array of values
  var idx = ctx.duk_push_array()
  for i, v in values:
    ctx.duk_push_number(v.duk_double_t)
    discard ctx.duk_put_prop_index(idx, i.duk_uarridx_t)
  discard ctx.duk_put_global_string(key)

proc `[]`*(ctx: DTContext, key:string): float {.inline.} =
    if ctx.duk_get_global_string(key) == 0:
        raise newException(KeyError, "couldn't find key:" & key)
    result = ctx.duk_get_number(-1).float

proc compile*(ctx: DTContext, expression: string): Dukexpr {.inline.} =
  ## compile an expression to be used later. This is about 10X faster than eval_string
  result = Dukexpr(ctx:ctx, expr: expression)
  if ctx.duk_pcompile_string(0, expression) != 0:
    var err = ctx.duk_safe_to_string(-1)
    raise newException(ValueError, $err)
  result.vptr = ctx.duk_get_heapptr(-1)

proc check*(d:Dukexpr): bool {.inline.} =
  ## evaluate a (previously compiled) boolean expression in the current context
  discard d.ctx.duk_push_heapptr(d.vptr)
  if d.ctx.duk_pcall(0) != 0:
    var err = d.ctx.duk_safe_to_string(-1)
    raise newException(ValueError, $err)
  result = d.ctx.duk_get_boolean(-1)
  d.ctx.pop()

proc check*(ctx: DTContext, expression: string): bool {.inline.} =
    ## evaluate the expression in the current context
    ctx.duk_eval_string(expression)
    result = ctx.duk_get_boolean(-1)
    ctx.pop()

proc newObject*(ctx:DTContext, name: string): Duko =
  ## create a new object.
  result = Duko(ctx: ctx, name: name)
  discard ctx.duk_push_object()
  result.vptr = ctx.duk_get_heapptr(-1)
  doAssert result.ctx.duk_put_global_string(name)

proc `[]=`*(o:Duko, key:string, value: SomeFloat) {.inline.} =
    ## set the property at key to a value
    var idx = o.ctx.duk_push_heapptr(o.vptr)
    o.ctx.duk_push_number(value.duk_double_t)
    doAssert o.ctx.duk_put_prop_string(idx, key)
    o.ctx.pop()

proc `[]=`*(o:Duko, key:string, value: SomeInteger) {.inline.} =
    ## set the property at key to a value
    var idx = o.ctx.duk_push_heapptr(o.vptr)
    o.ctx.duk_push_int(value.duk_int_t)
    if not o.ctx.duk_put_prop_string(idx, key):
      quit "problem setting:" & key & " -> " & $value
    o.ctx.pop()

proc alias*(o: Duko, copyname:string): Duko {.inline, discardable.} =
  ## create an alias of a Duko so it can be reference as another name in the javascript.
  doAssert o.ctx.duk_push_heapptr(o.vptr) >= 0
  result = Duko(ctx: o.ctx, name:copyname)
  result.vptr = o.vptr
  doAssert result.ctx.duk_put_global_string(copyname)

proc clear*(o: var Duko) {.inline.} =
  # TODO make this more efficient
  #o.ctx.duk_eval_string(o.name & "= null")
  discard o.ctx.duk_push_object()
  o.vptr = o.ctx.duk_get_heapptr(-1)
  doAssert o.ctx.duk_put_global_string(o.name)

proc `[]=`*(o: Duko, key: string, values: seq[SomeNumber]) {.inline.} =
  var idx = o.ctx.duk_push_heapptr(o.vptr)
  var arr_idx = o.ctx.duk_push_array()
  for i, v in values:
    o.ctx.duk_push_number(v.duk_double_t)
    discard o.ctx.duk_put_prop_index(arr_idx, i.duk_uarridx_t)
  doAssert o.ctx.duk_put_prop_string(idx, key)
  o.ctx.pop()

proc `[]`*(o: Duko, key:string): float {.inline.} =
    var idx = o.ctx.duk_push_heapptr(o.vptr)
    discard o.ctx.duk_push_string(key)
    doAssert o.ctx.duk_get_prop(idx)
    result = o.ctx.duk_get_number(-1).float
    o.ctx.duk_pop_n(2)

when isMainModule:
  import unittest
  var my_fatal: duk_fatal_function = (proc (udata: pointer, msg:cstring) {.stdcall.} =
      quit $msg
  )

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
      ctx.duk_destroy_heap();

    test "clear":
      var ctx = duk_create_heap_default()
      var obj = ctx.newObject("obj")
      obj["asdf"] = 1235
      obj["ddd"] = 22
      check obj["asdf"] == 1235
      obj.clear()
      ctx.duk_eval_string("obj.asdf")
      check ctx.duk_is_undefined(-1)

      ctx.duk_destroy_heap();

    test "speed":
      var ctx = duk_create_heap_default()
      var kid = ctx.newObject("kid")

      var t = cpuTime()
      when defined(release):
        var tries = 500_000
      else:
        var tries = 50_000

      var success = 0
      for i in 0..tries:
          kid["dp"] = i
          ctx["mom"] = i.float
          ctx["dad"] = 23
          ctx["proband"] = i.float
          kid["sdf"] = 22.2
          kid["xxx"] = 33.3 + i.float
          kid["yy"] = 33.4 + i.float

          if ctx.check("kid.dp < 500 && dad > 21 && mom == kid.dp"):
            success.inc
          kid.clear

      check success == 500
      echo (cpuTime() - t) / (tries / 1000000), " seconds per million evaluations"
      ctx.duk_destroy_heap();

    test "compiled speed":
      var ctx = duk_create_heap(nil, nil, nil, nil, my_fatal)
      var expr = "kid.dp < 500 && dad > 21" & " && mom == kid.dp"
      var e = ctx.compile(expr)
      var kid = ctx.newObject("kid")

      var t = cpuTime()
      when defined(release):
        var tries = 5_000_000
      else:
        var tries = 1_000_000

      var success = 0
      for i in 0..tries:
          kid["dp"] = i
          ctx["mom"] = i.float
          ctx["dad"] = 23
          ctx["proband"] = i.float
          kid["sdf"] = 22.2
          kid["xxx"] = 33.3 + i.float
          kid["yy"] = 33.4 + i.float
          var kid2 = kid.alias("kid2")
          check kid2.vptr == kid.vptr

          if e.check():
            success.inc
          kid.clear

      check success == 500
      echo (cpuTime() - t) / (tries / 1000000), " seconds per million evaluations"
      ctx.duk_destroy_heap();

    test "that alias works":
      var ctx = duk_create_heap_default()
      var kid = ctx.newObject("kid")
      for i in 0..100:
        kid["some" & $i] = i

      var kid2 = kid.alias("kid2")
      check kid2["some22"] == 22.0
      check kid["some22"] == 22.0
      kid2["some22"] = 44.0
      check kid["some22"] == 44.0
      ctx.duk_eval_string("kid.some22")
      check ctx.duk_get_number(-1) == 44.0
      ctx.duk_eval_string("kid2.some22")
      check ctx.duk_get_number(-1) == 44.0
      ctx.duk_destroy_heap();

    proc addo(obj: Duko) =
        for i in 0..200:
          obj["attr" & $i] = i.float
          obj["attrs" & $i] = @[i.float, i.float*2]

    test "many objects":
      var ctx = duk_create_heap(nil, nil, nil, nil, my_fatal)
      var objs = newSeq[Duko]()
      for i in 0..6000:
        objs.add(ctx.newObject("sample" & $i))
      for i in 0..10:
        for obj in objs.mitems:
          obj.clear
          obj.addo
        for obj in objs:
          obj.addo
          obj.alias("kid")
