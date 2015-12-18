var date_select = true;
!function(){function n(n,t){return t>n?-1:n>t?1:n>=t?0:0/0}function t(n){return null!=n&&!isNaN(n)}function e(n){return{left:function(t,e,r,u){for(arguments.length<3&&(r=0),arguments.length<4&&(u=t.length);u>r;){var i=r+u>>>1;n(t[i],e)<0?r=i+1:u=i}return r},right:function(t,e,r,u){for(arguments.length<3&&(r=0),arguments.length<4&&(u=t.length);u>r;){var i=r+u>>>1;n(t[i],e)>0?u=i:r=i+1}return r}}}function r(n){return n.length}function u(n){for(var t=1;n*t%1;)t*=10;return t}function i(n,t){try{for(var e in t)Object.defineProperty(n.prototype,e,{value:t[e],enumerable:!1})}catch(r){n.prototype=t}}function o(){}function a(n){return ia+n in this}function c(n){return n=ia+n,n in this&&delete this[n]}function s(){var n=[];return this.forEach(function(t){n.push(t)}),n}function l(){var n=0;for(var t in this)t.charCodeAt(0)===oa&&++n;return n}function f(){for(var n in this)if(n.charCodeAt(0)===oa)return!1;return!0}function h(){}function g(n,t,e){return function(){var r=e.apply(t,arguments);return r===t?n:r}}function p(n,t){if(t in n)return t;t=t.charAt(0).toUpperCase()+t.substring(1);for(var e=0,r=aa.length;r>e;++e){var u=aa[e]+t;if(u in n)return u}}function v(){}function d(){}function m(n){function t(){for(var t,r=e,u=-1,i=r.length;++u<i;)(t=r[u].on)&&t.apply(this,arguments);return n}var e=[],r=new o;return t.on=function(t,u){var i,o=r.get(t);return arguments.length<2?o&&o.on:(o&&(o.on=null,e=e.slice(0,i=e.indexOf(o)).concat(e.slice(i+1)),r.remove(t)),u&&e.push(r.set(t,{on:u})),n)},t}function y(){Zo.event.preventDefault()}function x(){for(var n,t=Zo.event;n=t.sourceEvent;)t=n;return t}function M(n){for(var t=new d,e=0,r=arguments.length;++e<r;)t[arguments[e]]=m(t);return t.of=function(e,r){return function(u){try{var i=u.sourceEvent=Zo.event;u.target=n,Zo.event=u,t[u.type].apply(e,r)}finally{Zo.event=i}}},t}function _(n){return sa(n,pa),n}function b(n){return"function"==typeof n?n:function(){return la(n,this)}}function w(n){return"function"==typeof n?n:function(){return fa(n,this)}}function S(n,t){function e(){this.removeAttribute(n)}function r(){this.removeAttributeNS(n.space,n.local)}function u(){this.setAttribute(n,t)}function i(){this.setAttributeNS(n.space,n.local,t)}function o(){var e=t.apply(this,arguments);null==e?this.removeAttribute(n):this.setAttribute(n,e)}function a(){var e=t.apply(this,arguments);null==e?this.removeAttributeNS(n.space,n.local):this.setAttributeNS(n.space,n.local,e)}return n=Zo.ns.qualify(n),null==t?n.local?r:e:"function"==typeof t?n.local?a:o:n.local?i:u}function k(n){return n.trim().replace(/\s+/g," ")}function E(n){return new RegExp("(?:^|\\s+)"+Zo.requote(n)+"(?:\\s+|$)","g")}function A(n){return(n+"").trim().split(/^|\s+/)}function C(n,t){function e(){for(var e=-1;++e<u;)n[e](this,t)}function r(){for(var e=-1,r=t.apply(this,arguments);++e<u;)n[e](this,r)}n=A(n).map(N);var u=n.length;return"function"==typeof t?r:e}function N(n){var t=E(n);return function(e,r){if(u=e.classList)return r?u.add(n):u.remove(n);var u=e.getAttribute("class")||"";r?(t.lastIndex=0,t.test(u)||e.setAttribute("class",k(u+" "+n))):e.setAttribute("class",k(u.replace(t," ")))}}function z(n,t,e){function r(){this.style.removeProperty(n)}function u(){this.style.setProperty(n,t,e)}function i(){var r=t.apply(this,arguments);null==r?this.style.removeProperty(n):this.style.setProperty(n,r,e)}return null==t?r:"function"==typeof t?i:u}function L(n,t){function e(){delete this[n]}function r(){this[n]=t}function u(){var e=t.apply(this,arguments);null==e?delete this[n]:this[n]=e}return null==t?e:"function"==typeof t?u:r}function T(n){return"function"==typeof n?n:(n=Zo.ns.qualify(n)).local?function(){return this.ownerDocument.createElementNS(n.space,n.local)}:function(){return this.ownerDocument.createElementNS(this.namespaceURI,n)}}function q(n){return{__data__:n}}function R(n){return function(){return ga(this,n)}}function D(t){return arguments.length||(t=n),function(n,e){return n&&e?t(n.__data__,e.__data__):!n-!e}}function P(n,t){for(var e=0,r=n.length;r>e;e++)for(var u,i=n[e],o=0,a=i.length;a>o;o++)(u=i[o])&&t(u,o,e);return n}function U(n){return sa(n,da),n}function j(n){var t,e;return function(r,u,i){var o,a=n[i].update,c=a.length;for(i!=e&&(e=i,t=0),u>=t&&(t=u+1);!(o=a[t])&&++t<c;);return o}}function H(){var n=this.__transition__;n&&++n.active}function F(n,t,e){function r(){var t=this[o];t&&(this.removeEventListener(n,t,t.$),delete this[o])}function u(){var u=c(t,Xo(arguments));r.call(this),this.addEventListener(n,this[o]=u,u.$=e),u._=t}function i(){var t,e=new RegExp("^__on([^.]+)"+Zo.requote(n)+"$");for(var r in this)if(t=r.match(e)){var u=this[r];this.removeEventListener(t[1],u,u.$),delete this[r]}}var o="__on"+n,a=n.indexOf("."),c=O;a>0&&(n=n.substring(0,a));var s=ya.get(n);return s&&(n=s,c=Y),a?t?u:r:t?v:i}function O(n,t){return function(e){var r=Zo.event;Zo.event=e,t[0]=this.__data__;try{n.apply(this,t)}finally{Zo.event=r}}}function Y(n,t){var e=O(n,t);return function(n){var t=this,r=n.relatedTarget;r&&(r===t||8&r.compareDocumentPosition(t))||e.call(t,n)}}function I(){var n=".dragsuppress-"+ ++Ma,t="click"+n,e=Zo.select(Wo).on("touchmove"+n,y).on("dragstart"+n,y).on("selectstart"+n,y);if(xa){var r=Bo.style,u=r[xa];r[xa]="none"}return function(i){function o(){e.on(t,null)}e.on(n,null),xa&&(r[xa]=u),i&&(e.on(t,function(){y(),o()},!0),setTimeout(o,0))}}function Z(n,t){t.changedTouches&&(t=t.changedTouches[0]);var e=n.ownerSVGElement||n;if(e.createSVGPoint){var r=e.createSVGPoint();if(0>_a&&(Wo.scrollX||Wo.scrollY)){e=Zo.select("body").append("svg").style({position:"absolute",top:0,left:0,margin:0,padding:0,border:"none"},"important");var u=e[0][0].getScreenCTM();_a=!(u.f||u.e),e.remove()}return _a?(r.x=t.pageX,r.y=t.pageY):(r.x=t.clientX,r.y=t.clientY),r=r.matrixTransform(n.getScreenCTM().inverse()),[r.x,r.y]}var i=n.getBoundingClientRect();return[t.clientX-i.left-n.clientLeft,t.clientY-i.top-n.clientTop]}function V(){return Zo.event.changedTouches[0].identifier}function X(){return Zo.event.target}function $(){return Wo}function B(n){return n>0?1:0>n?-1:0}function W(n,t,e){return(t[0]-n[0])*(e[1]-n[1])-(t[1]-n[1])*(e[0]-n[0])}function J(n){return n>1?0:-1>n?ba:Math.acos(n)}function G(n){return n>1?Sa:-1>n?-Sa:Math.asin(n)}function K(n){return((n=Math.exp(n))-1/n)/2}function Q(n){return((n=Math.exp(n))+1/n)/2}function nt(n){return((n=Math.exp(2*n))-1)/(n+1)}function tt(n){return(n=Math.sin(n/2))*n}function et(){}function rt(n,t,e){return this instanceof rt?(this.h=+n,this.s=+t,void(this.l=+e)):arguments.length<2?n instanceof rt?new rt(n.h,n.s,n.l):mt(""+n,yt,rt):new rt(n,t,e)}function ut(n,t,e){function r(n){return n>360?n-=360:0>n&&(n+=360),60>n?i+(o-i)*n/60:180>n?o:240>n?i+(o-i)*(240-n)/60:i}function u(n){return Math.round(255*r(n))}var i,o;return n=isNaN(n)?0:(n%=360)<0?n+360:n,t=isNaN(t)?0:0>t?0:t>1?1:t,e=0>e?0:e>1?1:e,o=.5>=e?e*(1+t):e+t-e*t,i=2*e-o,new gt(u(n+120),u(n),u(n-120))}function it(n,t,e){return this instanceof it?(this.h=+n,this.c=+t,void(this.l=+e)):arguments.length<2?n instanceof it?new it(n.h,n.c,n.l):n instanceof at?st(n.l,n.a,n.b):st((n=xt((n=Zo.rgb(n)).r,n.g,n.b)).l,n.a,n.b):new it(n,t,e)}function ot(n,t,e){return isNaN(n)&&(n=0),isNaN(t)&&(t=0),new at(e,Math.cos(n*=Aa)*t,Math.sin(n)*t)}function at(n,t,e){return this instanceof at?(this.l=+n,this.a=+t,void(this.b=+e)):arguments.length<2?n instanceof at?new at(n.l,n.a,n.b):n instanceof it?ot(n.l,n.c,n.h):xt((n=gt(n)).r,n.g,n.b):new at(n,t,e)}function ct(n,t,e){var r=(n+16)/116,u=r+t/500,i=r-e/200;return u=lt(u)*ja,r=lt(r)*Ha,i=lt(i)*Fa,new gt(ht(3.2404542*u-1.5371385*r-.4985314*i),ht(-.969266*u+1.8760108*r+.041556*i),ht(.0556434*u-.2040259*r+1.0572252*i))}function st(n,t,e){return n>0?new it(Math.atan2(e,t)*Ca,Math.sqrt(t*t+e*e),n):new it(0/0,0/0,n)}function lt(n){return n>.206893034?n*n*n:(n-4/29)/7.787037}function ft(n){return n>.008856?Math.pow(n,1/3):7.787037*n+4/29}function ht(n){return Math.round(255*(.00304>=n?12.92*n:1.055*Math.pow(n,1/2.4)-.055))}function gt(n,t,e){return this instanceof gt?(this.r=~~n,this.g=~~t,void(this.b=~~e)):arguments.length<2?n instanceof gt?new gt(n.r,n.g,n.b):mt(""+n,gt,ut):new gt(n,t,e)}function pt(n){return new gt(n>>16,255&n>>8,255&n)}function vt(n){return pt(n)+""}function dt(n){return 16>n?"0"+Math.max(0,n).toString(16):Math.min(255,n).toString(16)}function mt(n,t,e){var r,u,i,o=0,a=0,c=0;if(r=/([a-z]+)\((.*)\)/i.exec(n))switch(u=r[2].split(","),r[1]){case"hsl":return e(parseFloat(u[0]),parseFloat(u[1])/100,parseFloat(u[2])/100);case"rgb":return t(_t(u[0]),_t(u[1]),_t(u[2]))}return(i=Ia.get(n))?t(i.r,i.g,i.b):(null==n||"#"!==n.charAt(0)||isNaN(i=parseInt(n.substring(1),16))||(4===n.length?(o=(3840&i)>>4,o=o>>4|o,a=240&i,a=a>>4|a,c=15&i,c=c<<4|c):7===n.length&&(o=(16711680&i)>>16,a=(65280&i)>>8,c=255&i)),t(o,a,c))}function yt(n,t,e){var r,u,i=Math.min(n/=255,t/=255,e/=255),o=Math.max(n,t,e),a=o-i,c=(o+i)/2;return a?(u=.5>c?a/(o+i):a/(2-o-i),r=n==o?(t-e)/a+(e>t?6:0):t==o?(e-n)/a+2:(n-t)/a+4,r*=60):(r=0/0,u=c>0&&1>c?0:r),new rt(r,u,c)}function xt(n,t,e){n=Mt(n),t=Mt(t),e=Mt(e);var r=ft((.4124564*n+.3575761*t+.1804375*e)/ja),u=ft((.2126729*n+.7151522*t+.072175*e)/Ha),i=ft((.0193339*n+.119192*t+.9503041*e)/Fa);return at(116*u-16,500*(r-u),200*(u-i))}function Mt(n){return(n/=255)<=.04045?n/12.92:Math.pow((n+.055)/1.055,2.4)}function _t(n){var t=parseFloat(n);return"%"===n.charAt(n.length-1)?Math.round(2.55*t):t}function bt(n){return"function"==typeof n?n:function(){return n}}function wt(n){return n}function St(n){return function(t,e,r){return 2===arguments.length&&"function"==typeof e&&(r=e,e=null),kt(t,e,n,r)}}function kt(n,t,e,r){function u(){var n,t=c.status;if(!t&&c.responseText||t>=200&&300>t||304===t){try{n=e.call(i,c)}catch(r){return o.error.call(i,r),void 0}o.load.call(i,n)}else o.error.call(i,c)}var i={},o=Zo.dispatch("beforesend","progress","load","error"),a={},c=new XMLHttpRequest,s=null;return!Wo.XDomainRequest||"withCredentials"in c||!/^(http(s)?:)?\/\//.test(n)||(c=new XDomainRequest),"onload"in c?c.onload=c.onerror=u:c.onreadystatechange=function(){c.readyState>3&&u()},c.onprogress=function(n){var t=Zo.event;Zo.event=n;try{o.progress.call(i,c)}finally{Zo.event=t}},i.header=function(n,t){return n=(n+"").toLowerCase(),arguments.length<2?a[n]:(null==t?delete a[n]:a[n]=t+"",i)},i.mimeType=function(n){return arguments.length?(t=null==n?null:n+"",i):t},i.responseType=function(n){return arguments.length?(s=n,i):s},i.response=function(n){return e=n,i},["get","post"].forEach(function(n){i[n]=function(){return i.send.apply(i,[n].concat(Xo(arguments)))}}),i.send=function(e,r,u){if(2===arguments.length&&"function"==typeof r&&(u=r,r=null),c.open(e,n,!0),null==t||"accept"in a||(a.accept=t+",*/*"),c.setRequestHeader)for(var l in a)c.setRequestHeader(l,a[l]);return null!=t&&c.overrideMimeType&&c.overrideMimeType(t),null!=s&&(c.responseType=s),null!=u&&i.on("error",u).on("load",function(n){u(null,n)}),o.beforesend.call(i,c),c.send(null==r?null:r),i},i.abort=function(){return c.abort(),i},Zo.rebind(i,o,"on"),null==r?i:i.get(Et(r))}function Et(n){return 1===n.length?function(t,e){n(null==t?e:null)}:n}function At(){var n=Ct(),t=Nt()-n;t>24?(isFinite(t)&&(clearTimeout($a),$a=setTimeout(At,t)),Xa=0):(Xa=1,Wa(At))}function Ct(){var n=Date.now();for(Ba=Za;Ba;)n>=Ba.t&&(Ba.f=Ba.c(n-Ba.t)),Ba=Ba.n;return n}function Nt(){for(var n,t=Za,e=1/0;t;)t.f?t=n?n.n=t.n:Za=t.n:(t.t<e&&(e=t.t),t=(n=t).n);return Va=n,e}function zt(n,t){return t-(n?Math.ceil(Math.log(n)/Math.LN10):1)}function Lt(n,t){var e=Math.pow(10,3*ua(8-t));return{scale:t>8?function(n){return n/e}:function(n){return n*e},symbol:n}}function Tt(n){var t=n.decimal,e=n.thousands,r=n.grouping,u=n.currency,i=r?function(n){for(var t=n.length,u=[],i=0,o=r[0];t>0&&o>0;)u.push(n.substring(t-=o,t+o)),o=r[i=(i+1)%r.length];return u.reverse().join(e)}:wt;return function(n){var e=Ga.exec(n),r=e[1]||" ",o=e[2]||">",a=e[3]||"",c=e[4]||"",s=e[5],l=+e[6],f=e[7],h=e[8],g=e[9],p=1,v="",d="",m=!1;switch(h&&(h=+h.substring(1)),(s||"0"===r&&"="===o)&&(s=r="0",o="=",f&&(l-=Math.floor((l-1)/4))),g){case"n":f=!0,g="g";break;case"%":p=100,d="%",g="f";break;case"p":p=100,d="%",g="r";break;case"b":case"o":case"x":case"X":"#"===c&&(v="0"+g.toLowerCase());case"c":case"d":m=!0,h=0;break;case"s":p=-1,g="r"}"$"===c&&(v=u[0],d=u[1]),"r"!=g||h||(g="g"),null!=h&&("g"==g?h=Math.max(1,Math.min(21,h)):("e"==g||"f"==g)&&(h=Math.max(0,Math.min(20,h)))),g=Ka.get(g)||qt;var y=s&&f;return function(n){var e=d;if(m&&n%1)return"";var u=0>n||0===n&&0>1/n?(n=-n,"-"):a;if(0>p){var c=Zo.formatPrefix(n,h);n=c.scale(n),e=c.symbol+d}else n*=p;n=g(n,h);var x=n.lastIndexOf("."),M=0>x?n:n.substring(0,x),_=0>x?"":t+n.substring(x+1);!s&&f&&(M=i(M));var b=v.length+M.length+_.length+(y?0:u.length),w=l>b?new Array(b=l-b+1).join(r):"";return y&&(M=i(w+M)),u+=v,n=M+_,("<"===o?u+n+w:">"===o?w+u+n:"^"===o?w.substring(0,b>>=1)+u+n+w.substring(b):u+(y?n:w+n))+e}}}function qt(n){return n+""}function Rt(){this._=new Date(arguments.length>1?Date.UTC.apply(this,arguments):arguments[0])}function Dt(n,t,e){function r(t){var e=n(t),r=i(e,1);return r-t>t-e?e:r}function u(e){return t(e=n(new nc(e-1)),1),e}function i(n,e){return t(n=new nc(+n),e),n}function o(n,r,i){var o=u(n),a=[];if(i>1)for(;r>o;)e(o)%i||a.push(new Date(+o)),t(o,1);else for(;r>o;)a.push(new Date(+o)),t(o,1);return a}function a(n,t,e){try{nc=Rt;var r=new Rt;return r._=n,o(r,t,e)}finally{nc=Date}}n.floor=n,n.round=r,n.ceil=u,n.offset=i,n.range=o;var c=n.utc=Pt(n);return c.floor=c,c.round=Pt(r),c.ceil=Pt(u),c.offset=Pt(i),c.range=a,n}function Pt(n){return function(t,e){try{nc=Rt;var r=new Rt;return r._=t,n(r,e)._}finally{nc=Date}}}function Ut(n){function t(n){function t(t){for(var e,u,i,o=[],a=-1,c=0;++a<r;)37===n.charCodeAt(a)&&(o.push(n.substring(c,a)),null!=(u=ec[e=n.charAt(++a)])&&(e=n.charAt(++a)),(i=C[e])&&(e=i(t,null==u?"e"===e?" ":"0":u)),o.push(e),c=a+1);return o.push(n.substring(c,a)),o.join("")}var r=n.length;return t.parse=function(t){var r={y:1900,m:0,d:1,H:0,M:0,S:0,L:0,Z:null},u=e(r,n,t,0);if(u!=t.length)return null;"p"in r&&(r.H=r.H%12+12*r.p);var i=null!=r.Z&&nc!==Rt,o=new(i?Rt:nc);return"j"in r?o.setFullYear(r.y,0,r.j):"w"in r&&("W"in r||"U"in r)?(o.setFullYear(r.y,0,1),o.setFullYear(r.y,0,"W"in r?(r.w+6)%7+7*r.W-(o.getDay()+5)%7:r.w+7*r.U-(o.getDay()+6)%7)):o.setFullYear(r.y,r.m,r.d),o.setHours(r.H+Math.floor(r.Z/100),r.M+r.Z%100,r.S,r.L),i?o._:o},t.toString=function(){return n},t}function e(n,t,e,r){for(var u,i,o,a=0,c=t.length,s=e.length;c>a;){if(r>=s)return-1;if(u=t.charCodeAt(a++),37===u){if(o=t.charAt(a++),i=N[o in ec?t.charAt(a++):o],!i||(r=i(n,e,r))<0)return-1}else if(u!=e.charCodeAt(r++))return-1}return r}function r(n,t,e){b.lastIndex=0;var r=b.exec(t.substring(e));return r?(n.w=w.get(r[0].toLowerCase()),e+r[0].length):-1}function u(n,t,e){M.lastIndex=0;var r=M.exec(t.substring(e));return r?(n.w=_.get(r[0].toLowerCase()),e+r[0].length):-1}function i(n,t,e){E.lastIndex=0;var r=E.exec(t.substring(e));return r?(n.m=A.get(r[0].toLowerCase()),e+r[0].length):-1}function o(n,t,e){S.lastIndex=0;var r=S.exec(t.substring(e));return r?(n.m=k.get(r[0].toLowerCase()),e+r[0].length):-1}function a(n,t,r){return e(n,C.c.toString(),t,r)}function c(n,t,r){return e(n,C.x.toString(),t,r)}function s(n,t,r){return e(n,C.X.toString(),t,r)}function l(n,t,e){var r=x.get(t.substring(e,e+=2).toLowerCase());return null==r?-1:(n.p=r,e)}var f=n.dateTime,h=n.date,g=n.time,p=n.periods,v=n.days,d=n.shortDays,m=n.months,y=n.shortMonths;t.utc=function(n){function e(n){try{nc=Rt;var t=new nc;return t._=n,r(t)}finally{nc=Date}}var r=t(n);return e.parse=function(n){try{nc=Rt;var t=r.parse(n);return t&&t._}finally{nc=Date}},e.toString=r.toString,e},t.multi=t.utc.multi=re;var x=Zo.map(),M=Ht(v),_=Ft(v),b=Ht(d),w=Ft(d),S=Ht(m),k=Ft(m),E=Ht(y),A=Ft(y);p.forEach(function(n,t){x.set(n.toLowerCase(),t)});var C={a:function(n){return d[n.getDay()]},A:function(n){return v[n.getDay()]},b:function(n){return y[n.getMonth()]},B:function(n){return m[n.getMonth()]},c:t(f),d:function(n,t){return jt(n.getDate(),t,2)},e:function(n,t){return jt(n.getDate(),t,2)},H:function(n,t){return jt(n.getHours(),t,2)},I:function(n,t){return jt(n.getHours()%12||12,t,2)},j:function(n,t){return jt(1+Qa.dayOfYear(n),t,3)},L:function(n,t){return jt(n.getMilliseconds(),t,3)},m:function(n,t){return jt(n.getMonth()+1,t,2)},M:function(n,t){return jt(n.getMinutes(),t,2)},p:function(n){return p[+(n.getHours()>=12)]},S:function(n,t){return jt(n.getSeconds(),t,2)},U:function(n,t){return jt(Qa.sundayOfYear(n),t,2)},w:function(n){return n.getDay()},W:function(n,t){return jt(Qa.mondayOfYear(n),t,2)},x:t(h),X:t(g),y:function(n,t){return jt(n.getFullYear()%100,t,2)},Y:function(n,t){return jt(n.getFullYear()%1e4,t,4)},Z:te,"%":function(){return"%"}},N={a:r,A:u,b:i,B:o,c:a,d:Wt,e:Wt,H:Gt,I:Gt,j:Jt,L:ne,m:Bt,M:Kt,p:l,S:Qt,U:Yt,w:Ot,W:It,x:c,X:s,y:Vt,Y:Zt,Z:Xt,"%":ee};return t}function jt(n,t,e){var r=0>n?"-":"",u=(r?-n:n)+"",i=u.length;return r+(e>i?new Array(e-i+1).join(t)+u:u)}function Ht(n){return new RegExp("^(?:"+n.map(Zo.requote).join("|")+")","i")}function Ft(n){for(var t=new o,e=-1,r=n.length;++e<r;)t.set(n[e].toLowerCase(),e);return t}function Ot(n,t,e){rc.lastIndex=0;var r=rc.exec(t.substring(e,e+1));return r?(n.w=+r[0],e+r[0].length):-1}function Yt(n,t,e){rc.lastIndex=0;var r=rc.exec(t.substring(e));return r?(n.U=+r[0],e+r[0].length):-1}function It(n,t,e){rc.lastIndex=0;var r=rc.exec(t.substring(e));return r?(n.W=+r[0],e+r[0].length):-1}function Zt(n,t,e){rc.lastIndex=0;var r=rc.exec(t.substring(e,e+4));return r?(n.y=+r[0],e+r[0].length):-1}function Vt(n,t,e){rc.lastIndex=0;var r=rc.exec(t.substring(e,e+2));return r?(n.y=$t(+r[0]),e+r[0].length):-1}function Xt(n,t,e){return/^[+-]\d{4}$/.test(t=t.substring(e,e+5))?(n.Z=-t,e+5):-1}function $t(n){return n+(n>68?1900:2e3)}function Bt(n,t,e){rc.lastIndex=0;var r=rc.exec(t.substring(e,e+2));return r?(n.m=r[0]-1,e+r[0].length):-1}function Wt(n,t,e){rc.lastIndex=0;var r=rc.exec(t.substring(e,e+2));return r?(n.d=+r[0],e+r[0].length):-1}function Jt(n,t,e){rc.lastIndex=0;var r=rc.exec(t.substring(e,e+3));return r?(n.j=+r[0],e+r[0].length):-1}function Gt(n,t,e){rc.lastIndex=0;var r=rc.exec(t.substring(e,e+2));return r?(n.H=+r[0],e+r[0].length):-1}function Kt(n,t,e){rc.lastIndex=0;var r=rc.exec(t.substring(e,e+2));return r?(n.M=+r[0],e+r[0].length):-1}function Qt(n,t,e){rc.lastIndex=0;var r=rc.exec(t.substring(e,e+2));return r?(n.S=+r[0],e+r[0].length):-1}function ne(n,t,e){rc.lastIndex=0;var r=rc.exec(t.substring(e,e+3));return r?(n.L=+r[0],e+r[0].length):-1}function te(n){var t=n.getTimezoneOffset(),e=t>0?"-":"+",r=~~(ua(t)/60),u=ua(t)%60;return e+jt(r,"0",2)+jt(u,"0",2)}function ee(n,t,e){uc.lastIndex=0;var r=uc.exec(t.substring(e,e+1));return r?e+r[0].length:-1}function re(n){for(var t=n.length,e=-1;++e<t;)n[e][0]=this(n[e][0]);return function(t){for(var e=0,r=n[e];!r[1](t);)r=n[++e];return r[0](t)}}function ue(){}function ie(n,t,e){var r=e.s=n+t,u=r-n,i=r-u;e.t=n-i+(t-u)}function oe(n,t){n&&cc.hasOwnProperty(n.type)&&cc[n.type](n,t)}function ae(n,t,e){var r,u=-1,i=n.length-e;for(t.lineStart();++u<i;)r=n[u],t.point(r[0],r[1],r[2]);t.lineEnd()}function ce(n,t){var e=-1,r=n.length;for(t.polygonStart();++e<r;)ae(n[e],t,1);t.polygonEnd()}function se(){function n(n,t){n*=Aa,t=t*Aa/2+ba/4;var e=n-r,o=e>=0?1:-1,a=o*e,c=Math.cos(t),s=Math.sin(t),l=i*s,f=u*c+l*Math.cos(a),h=l*o*Math.sin(a);lc.add(Math.atan2(h,f)),r=n,u=c,i=s}var t,e,r,u,i;fc.point=function(o,a){fc.point=n,r=(t=o)*Aa,u=Math.cos(a=(e=a)*Aa/2+ba/4),i=Math.sin(a)},fc.lineEnd=function(){n(t,e)}}function le(n){var t=n[0],e=n[1],r=Math.cos(e);return[r*Math.cos(t),r*Math.sin(t),Math.sin(e)]}function fe(n,t){return n[0]*t[0]+n[1]*t[1]+n[2]*t[2]}function he(n,t){return[n[1]*t[2]-n[2]*t[1],n[2]*t[0]-n[0]*t[2],n[0]*t[1]-n[1]*t[0]]}function ge(n,t){n[0]+=t[0],n[1]+=t[1],n[2]+=t[2]}function pe(n,t){return[n[0]*t,n[1]*t,n[2]*t]}function ve(n){var t=Math.sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);n[0]/=t,n[1]/=t,n[2]/=t}function de(n){return[Math.atan2(n[1],n[0]),G(n[2])]}function me(n,t){return ua(n[0]-t[0])<ka&&ua(n[1]-t[1])<ka}function ye(n,t){n*=Aa;var e=Math.cos(t*=Aa);xe(e*Math.cos(n),e*Math.sin(n),Math.sin(t))}function xe(n,t,e){++hc,pc+=(n-pc)/hc,vc+=(t-vc)/hc,dc+=(e-dc)/hc}function Me(){function n(n,u){n*=Aa;var i=Math.cos(u*=Aa),o=i*Math.cos(n),a=i*Math.sin(n),c=Math.sin(u),s=Math.atan2(Math.sqrt((s=e*c-r*a)*s+(s=r*o-t*c)*s+(s=t*a-e*o)*s),t*o+e*a+r*c);gc+=s,mc+=s*(t+(t=o)),yc+=s*(e+(e=a)),xc+=s*(r+(r=c)),xe(t,e,r)}var t,e,r;wc.point=function(u,i){u*=Aa;var o=Math.cos(i*=Aa);t=o*Math.cos(u),e=o*Math.sin(u),r=Math.sin(i),wc.point=n,xe(t,e,r)}}function _e(){wc.point=ye}function be(){function n(n,t){n*=Aa;var e=Math.cos(t*=Aa),o=e*Math.cos(n),a=e*Math.sin(n),c=Math.sin(t),s=u*c-i*a,l=i*o-r*c,f=r*a-u*o,h=Math.sqrt(s*s+l*l+f*f),g=r*o+u*a+i*c,p=h&&-J(g)/h,v=Math.atan2(h,g);Mc+=p*s,_c+=p*l,bc+=p*f,gc+=v,mc+=v*(r+(r=o)),yc+=v*(u+(u=a)),xc+=v*(i+(i=c)),xe(r,u,i)}var t,e,r,u,i;wc.point=function(o,a){t=o,e=a,wc.point=n,o*=Aa;var c=Math.cos(a*=Aa);r=c*Math.cos(o),u=c*Math.sin(o),i=Math.sin(a),xe(r,u,i)},wc.lineEnd=function(){n(t,e),wc.lineEnd=_e,wc.point=ye}}function we(){return!0}function Se(n,t,e,r,u){var i=[],o=[];if(n.forEach(function(n){if(!((t=n.length-1)<=0)){var t,e=n[0],r=n[t];if(me(e,r)){u.lineStart();for(var a=0;t>a;++a)u.point((e=n[a])[0],e[1]);return u.lineEnd(),void 0}var c=new Ee(e,n,null,!0),s=new Ee(e,null,c,!1);c.o=s,i.push(c),o.push(s),c=new Ee(r,n,null,!1),s=new Ee(r,null,c,!0),c.o=s,i.push(c),o.push(s)}}),o.sort(t),ke(i),ke(o),i.length){for(var a=0,c=e,s=o.length;s>a;++a)o[a].e=c=!c;for(var l,f,h=i[0];;){for(var g=h,p=!0;g.v;)if((g=g.n)===h)return;l=g.z,u.lineStart();do{if(g.v=g.o.v=!0,g.e){if(p)for(var a=0,s=l.length;s>a;++a)u.point((f=l[a])[0],f[1]);else r(g.x,g.n.x,1,u);g=g.n}else{if(p){l=g.p.z;for(var a=l.length-1;a>=0;--a)u.point((f=l[a])[0],f[1])}else r(g.x,g.p.x,-1,u);g=g.p}g=g.o,l=g.z,p=!p}while(!g.v);u.lineEnd()}}}function ke(n){if(t=n.length){for(var t,e,r=0,u=n[0];++r<t;)u.n=e=n[r],e.p=u,u=e;u.n=e=n[0],e.p=u}}function Ee(n,t,e,r){this.x=n,this.z=t,this.o=e,this.e=r,this.v=!1,this.n=this.p=null}function Ae(n,t,e,r){return function(u,i){function o(t,e){var r=u(t,e);n(t=r[0],e=r[1])&&i.point(t,e)}function a(n,t){var e=u(n,t);d.point(e[0],e[1])}function c(){y.point=a,d.lineStart()}function s(){y.point=o,d.lineEnd()}function l(n,t){v.push([n,t]);var e=u(n,t);M.point(e[0],e[1])}function f(){M.lineStart(),v=[]}function h(){l(v[0][0],v[0][1]),M.lineEnd();var n,t=M.clean(),e=x.buffer(),r=e.length;if(v.pop(),p.push(v),v=null,r)if(1&t){n=e[0];var u,r=n.length-1,o=-1;if(r>0){for(_||(i.polygonStart(),_=!0),i.lineStart();++o<r;)i.point((u=n[o])[0],u[1]);i.lineEnd()}}else r>1&&2&t&&e.push(e.pop().concat(e.shift())),g.push(e.filter(Ce))}var g,p,v,d=t(i),m=u.invert(r[0],r[1]),y={point:o,lineStart:c,lineEnd:s,polygonStart:function(){y.point=l,y.lineStart=f,y.lineEnd=h,g=[],p=[]},polygonEnd:function(){y.point=o,y.lineStart=c,y.lineEnd=s,g=Zo.merge(g);var n=Le(m,p);g.length?(_||(i.polygonStart(),_=!0),Se(g,ze,n,e,i)):n&&(_||(i.polygonStart(),_=!0),i.lineStart(),e(null,null,1,i),i.lineEnd()),_&&(i.polygonEnd(),_=!1),g=p=null},sphere:function(){i.polygonStart(),i.lineStart(),e(null,null,1,i),i.lineEnd(),i.polygonEnd()}},x=Ne(),M=t(x),_=!1;return y}}function Ce(n){return n.length>1}function Ne(){var n,t=[];return{lineStart:function(){t.push(n=[])},point:function(t,e){n.push([t,e])},lineEnd:v,buffer:function(){var e=t;return t=[],n=null,e},rejoin:function(){t.length>1&&t.push(t.pop().concat(t.shift()))}}}function ze(n,t){return((n=n.x)[0]<0?n[1]-Sa-ka:Sa-n[1])-((t=t.x)[0]<0?t[1]-Sa-ka:Sa-t[1])}function Le(n,t){var e=n[0],r=n[1],u=[Math.sin(e),-Math.cos(e),0],i=0,o=0;lc.reset();for(var a=0,c=t.length;c>a;++a){var s=t[a],l=s.length;if(l)for(var f=s[0],h=f[0],g=f[1]/2+ba/4,p=Math.sin(g),v=Math.cos(g),d=1;;){d===l&&(d=0),n=s[d];var m=n[0],y=n[1]/2+ba/4,x=Math.sin(y),M=Math.cos(y),_=m-h,b=_>=0?1:-1,w=b*_,S=w>ba,k=p*x;if(lc.add(Math.atan2(k*b*Math.sin(w),v*M+k*Math.cos(w))),i+=S?_+b*wa:_,S^h>=e^m>=e){var E=he(le(f),le(n));ve(E);var A=he(u,E);ve(A);var C=(S^_>=0?-1:1)*G(A[2]);(r>C||r===C&&(E[0]||E[1]))&&(o+=S^_>=0?1:-1)}if(!d++)break;h=m,p=x,v=M,f=n}}return(-ka>i||ka>i&&0>lc)^1&o}function Te(n){var t,e=0/0,r=0/0,u=0/0;return{lineStart:function(){n.lineStart(),t=1},point:function(i,o){var a=i>0?ba:-ba,c=ua(i-e);ua(c-ba)<ka?(n.point(e,r=(r+o)/2>0?Sa:-Sa),n.point(u,r),n.lineEnd(),n.lineStart(),n.point(a,r),n.point(i,r),t=0):u!==a&&c>=ba&&(ua(e-u)<ka&&(e-=u*ka),ua(i-a)<ka&&(i-=a*ka),r=qe(e,r,i,o),n.point(u,r),n.lineEnd(),n.lineStart(),n.point(a,r),t=0),n.point(e=i,r=o),u=a},lineEnd:function(){n.lineEnd(),e=r=0/0},clean:function(){return 2-t}}}function qe(n,t,e,r){var u,i,o=Math.sin(n-e);return ua(o)>ka?Math.atan((Math.sin(t)*(i=Math.cos(r))*Math.sin(e)-Math.sin(r)*(u=Math.cos(t))*Math.sin(n))/(u*i*o)):(t+r)/2}function Re(n,t,e,r){var u;if(null==n)u=e*Sa,r.point(-ba,u),r.point(0,u),r.point(ba,u),r.point(ba,0),r.point(ba,-u),r.point(0,-u),r.point(-ba,-u),r.point(-ba,0),r.point(-ba,u);else if(ua(n[0]-t[0])>ka){var i=n[0]<t[0]?ba:-ba;u=e*i/2,r.point(-i,u),r.point(0,u),r.point(i,u)}else r.point(t[0],t[1])}function De(n){function t(n,t){return Math.cos(n)*Math.cos(t)>i}function e(n){var e,i,c,s,l;return{lineStart:function(){s=c=!1,l=1},point:function(f,h){var g,p=[f,h],v=t(f,h),d=o?v?0:u(f,h):v?u(f+(0>f?ba:-ba),h):0;if(!e&&(s=c=v)&&n.lineStart(),v!==c&&(g=r(e,p),(me(e,g)||me(p,g))&&(p[0]+=ka,p[1]+=ka,v=t(p[0],p[1]))),v!==c)l=0,v?(n.lineStart(),g=r(p,e),n.point(g[0],g[1])):(g=r(e,p),n.point(g[0],g[1]),n.lineEnd()),e=g;else if(a&&e&&o^v){var m;d&i||!(m=r(p,e,!0))||(l=0,o?(n.lineStart(),n.point(m[0][0],m[0][1]),n.point(m[1][0],m[1][1]),n.lineEnd()):(n.point(m[1][0],m[1][1]),n.lineEnd(),n.lineStart(),n.point(m[0][0],m[0][1])))}!v||e&&me(e,p)||n.point(p[0],p[1]),e=p,c=v,i=d},lineEnd:function(){c&&n.lineEnd(),e=null},clean:function(){return l|(s&&c)<<1}}}function r(n,t,e){var r=le(n),u=le(t),o=[1,0,0],a=he(r,u),c=fe(a,a),s=a[0],l=c-s*s;if(!l)return!e&&n;var f=i*c/l,h=-i*s/l,g=he(o,a),p=pe(o,f),v=pe(a,h);ge(p,v);var d=g,m=fe(p,d),y=fe(d,d),x=m*m-y*(fe(p,p)-1);if(!(0>x)){var M=Math.sqrt(x),_=pe(d,(-m-M)/y);if(ge(_,p),_=de(_),!e)return _;var b,w=n[0],S=t[0],k=n[1],E=t[1];w>S&&(b=w,w=S,S=b);var A=S-w,C=ua(A-ba)<ka,N=C||ka>A;if(!C&&k>E&&(b=k,k=E,E=b),N?C?k+E>0^_[1]<(ua(_[0]-w)<ka?k:E):k<=_[1]&&_[1]<=E:A>ba^(w<=_[0]&&_[0]<=S)){var z=pe(d,(-m+M)/y);return ge(z,p),[_,de(z)]}}}function u(t,e){var r=o?n:ba-n,u=0;return-r>t?u|=1:t>r&&(u|=2),-r>e?u|=4:e>r&&(u|=8),u}var i=Math.cos(n),o=i>0,a=ua(i)>ka,c=sr(n,6*Aa);return Ae(t,e,c,o?[0,-n]:[-ba,n-ba])}function Pe(n,t,e,r){return function(u){var i,o=u.a,a=u.b,c=o.x,s=o.y,l=a.x,f=a.y,h=0,g=1,p=l-c,v=f-s;if(i=n-c,p||!(i>0)){if(i/=p,0>p){if(h>i)return;g>i&&(g=i)}else if(p>0){if(i>g)return;i>h&&(h=i)}if(i=e-c,p||!(0>i)){if(i/=p,0>p){if(i>g)return;i>h&&(h=i)}else if(p>0){if(h>i)return;g>i&&(g=i)}if(i=t-s,v||!(i>0)){if(i/=v,0>v){if(h>i)return;g>i&&(g=i)}else if(v>0){if(i>g)return;i>h&&(h=i)}if(i=r-s,v||!(0>i)){if(i/=v,0>v){if(i>g)return;i>h&&(h=i)}else if(v>0){if(h>i)return;g>i&&(g=i)}return h>0&&(u.a={x:c+h*p,y:s+h*v}),1>g&&(u.b={x:c+g*p,y:s+g*v}),u}}}}}}function Ue(n,t,e,r){function u(r,u){return ua(r[0]-n)<ka?u>0?0:3:ua(r[0]-e)<ka?u>0?2:1:ua(r[1]-t)<ka?u>0?1:0:u>0?3:2}function i(n,t){return o(n.x,t.x)}function o(n,t){var e=u(n,1),r=u(t,1);return e!==r?e-r:0===e?t[1]-n[1]:1===e?n[0]-t[0]:2===e?n[1]-t[1]:t[0]-n[0]}return function(a){function c(n){for(var t=0,e=d.length,r=n[1],u=0;e>u;++u)for(var i,o=1,a=d[u],c=a.length,s=a[0];c>o;++o)i=a[o],s[1]<=r?i[1]>r&&W(s,i,n)>0&&++t:i[1]<=r&&W(s,i,n)<0&&--t,s=i;return 0!==t}function s(i,a,c,s){var l=0,f=0;if(null==i||(l=u(i,c))!==(f=u(a,c))||o(i,a)<0^c>0){do s.point(0===l||3===l?n:e,l>1?r:t);while((l=(l+c+4)%4)!==f)}else s.point(a[0],a[1])}function l(u,i){return u>=n&&e>=u&&i>=t&&r>=i}function f(n,t){l(n,t)&&a.point(n,t)}function h(){N.point=p,d&&d.push(m=[]),S=!0,w=!1,_=b=0/0}function g(){v&&(p(y,x),M&&w&&A.rejoin(),v.push(A.buffer())),N.point=f,w&&a.lineEnd()}function p(n,t){n=Math.max(-kc,Math.min(kc,n)),t=Math.max(-kc,Math.min(kc,t));var e=l(n,t);if(d&&m.push([n,t]),S)y=n,x=t,M=e,S=!1,e&&(a.lineStart(),a.point(n,t));else if(e&&w)a.point(n,t);else{var r={a:{x:_,y:b},b:{x:n,y:t}};C(r)?(w||(a.lineStart(),a.point(r.a.x,r.a.y)),a.point(r.b.x,r.b.y),e||a.lineEnd(),k=!1):e&&(a.lineStart(),a.point(n,t),k=!1)}_=n,b=t,w=e}var v,d,m,y,x,M,_,b,w,S,k,E=a,A=Ne(),C=Pe(n,t,e,r),N={point:f,lineStart:h,lineEnd:g,polygonStart:function(){a=A,v=[],d=[],k=!0},polygonEnd:function(){a=E,v=Zo.merge(v);var t=c([n,r]),e=k&&t,u=v.length;(e||u)&&(a.polygonStart(),e&&(a.lineStart(),s(null,null,1,a),a.lineEnd()),u&&Se(v,i,t,s,a),a.polygonEnd()),v=d=m=null}};return N}}function je(n,t){function e(e,r){return e=n(e,r),t(e[0],e[1])}return n.invert&&t.invert&&(e.invert=function(e,r){return e=t.invert(e,r),e&&n.invert(e[0],e[1])}),e}function He(n){var t=0,e=ba/3,r=tr(n),u=r(t,e);return u.parallels=function(n){return arguments.length?r(t=n[0]*ba/180,e=n[1]*ba/180):[180*(t/ba),180*(e/ba)]},u}function Fe(n,t){function e(n,t){var e=Math.sqrt(i-2*u*Math.sin(t))/u;return[e*Math.sin(n*=u),o-e*Math.cos(n)]}var r=Math.sin(n),u=(r+Math.sin(t))/2,i=1+r*(2*u-r),o=Math.sqrt(i)/u;return e.invert=function(n,t){var e=o-t;return[Math.atan2(n,e)/u,G((i-(n*n+e*e)*u*u)/(2*u))]},e}function Oe(){function n(n,t){Ac+=u*n-r*t,r=n,u=t}var t,e,r,u;Tc.point=function(i,o){Tc.point=n,t=r=i,e=u=o},Tc.lineEnd=function(){n(t,e)}}function Ye(n,t){Cc>n&&(Cc=n),n>zc&&(zc=n),Nc>t&&(Nc=t),t>Lc&&(Lc=t)}function Ie(){function n(n,t){o.push("M",n,",",t,i)}function t(n,t){o.push("M",n,",",t),a.point=e}function e(n,t){o.push("L",n,",",t)}function r(){a.point=n}function u(){o.push("Z")}var i=Ze(4.5),o=[],a={point:n,lineStart:function(){a.point=t},lineEnd:r,polygonStart:function(){a.lineEnd=u},polygonEnd:function(){a.lineEnd=r,a.point=n},pointRadius:function(n){return i=Ze(n),a},result:function(){if(o.length){var n=o.join("");return o=[],n}}};return a}function Ze(n){return"m0,"+n+"a"+n+","+n+" 0 1,1 0,"+-2*n+"a"+n+","+n+" 0 1,1 0,"+2*n+"z"}function Ve(n,t){pc+=n,vc+=t,++dc}function Xe(){function n(n,r){var u=n-t,i=r-e,o=Math.sqrt(u*u+i*i);mc+=o*(t+n)/2,yc+=o*(e+r)/2,xc+=o,Ve(t=n,e=r)}var t,e;Rc.point=function(r,u){Rc.point=n,Ve(t=r,e=u)}}function $e(){Rc.point=Ve}function Be(){function n(n,t){var e=n-r,i=t-u,o=Math.sqrt(e*e+i*i);mc+=o*(r+n)/2,yc+=o*(u+t)/2,xc+=o,o=u*n-r*t,Mc+=o*(r+n),_c+=o*(u+t),bc+=3*o,Ve(r=n,u=t)}var t,e,r,u;Rc.point=function(i,o){Rc.point=n,Ve(t=r=i,e=u=o)},Rc.lineEnd=function(){n(t,e)}}function We(n){function t(t,e){n.moveTo(t,e),n.arc(t,e,o,0,wa)}function e(t,e){n.moveTo(t,e),a.point=r}function r(t,e){n.lineTo(t,e)}function u(){a.point=t}function i(){n.closePath()}var o=4.5,a={point:t,lineStart:function(){a.point=e},lineEnd:u,polygonStart:function(){a.lineEnd=i},polygonEnd:function(){a.lineEnd=u,a.point=t},pointRadius:function(n){return o=n,a},result:v};return a}function Je(n){function t(n){return(a?r:e)(n)}function e(t){return Qe(t,function(e,r){e=n(e,r),t.point(e[0],e[1])})}function r(t){function e(e,r){e=n(e,r),t.point(e[0],e[1])}function r(){x=0/0,S.point=i,t.lineStart()}function i(e,r){var i=le([e,r]),o=n(e,r);u(x,M,y,_,b,w,x=o[0],M=o[1],y=e,_=i[0],b=i[1],w=i[2],a,t),t.point(x,M)}function o(){S.point=e,t.lineEnd()}function c(){r(),S.point=s,S.lineEnd=l}function s(n,t){i(f=n,h=t),g=x,p=M,v=_,d=b,m=w,S.point=i}function l(){u(x,M,y,_,b,w,g,p,f,v,d,m,a,t),S.lineEnd=o,o()}var f,h,g,p,v,d,m,y,x,M,_,b,w,S={point:e,lineStart:r,lineEnd:o,polygonStart:function(){t.polygonStart(),S.lineStart=c},polygonEnd:function(){t.polygonEnd(),S.lineStart=r}};return S}function u(t,e,r,a,c,s,l,f,h,g,p,v,d,m){var y=l-t,x=f-e,M=y*y+x*x;if(M>4*i&&d--){var _=a+g,b=c+p,w=s+v,S=Math.sqrt(_*_+b*b+w*w),k=Math.asin(w/=S),E=ua(ua(w)-1)<ka||ua(r-h)<ka?(r+h)/2:Math.atan2(b,_),A=n(E,k),C=A[0],N=A[1],z=C-t,L=N-e,T=x*z-y*L;(T*T/M>i||ua((y*z+x*L)/M-.5)>.3||o>a*g+c*p+s*v)&&(u(t,e,r,a,c,s,C,N,E,_/=S,b/=S,w,d,m),m.point(C,N),u(C,N,E,_,b,w,l,f,h,g,p,v,d,m))}}var i=.5,o=Math.cos(30*Aa),a=16;
return t.precision=function(n){return arguments.length?(a=(i=n*n)>0&&16,t):Math.sqrt(i)},t}function Ge(n){var t=Je(function(t,e){return n([t*Ca,e*Ca])});return function(n){return er(t(n))}}function Ke(n){this.stream=n}function Qe(n,t){return{point:t,sphere:function(){n.sphere()},lineStart:function(){n.lineStart()},lineEnd:function(){n.lineEnd()},polygonStart:function(){n.polygonStart()},polygonEnd:function(){n.polygonEnd()}}}function nr(n){return tr(function(){return n})()}function tr(n){function t(n){return n=a(n[0]*Aa,n[1]*Aa),[n[0]*h+c,s-n[1]*h]}function e(n){return n=a.invert((n[0]-c)/h,(s-n[1])/h),n&&[n[0]*Ca,n[1]*Ca]}function r(){a=je(o=ir(m,y,x),i);var n=i(v,d);return c=g-n[0]*h,s=p+n[1]*h,u()}function u(){return l&&(l.valid=!1,l=null),t}var i,o,a,c,s,l,f=Je(function(n,t){return n=i(n,t),[n[0]*h+c,s-n[1]*h]}),h=150,g=480,p=250,v=0,d=0,m=0,y=0,x=0,M=Sc,_=wt,b=null,w=null;return t.stream=function(n){return l&&(l.valid=!1),l=er(M(o,f(_(n)))),l.valid=!0,l},t.clipAngle=function(n){return arguments.length?(M=null==n?(b=n,Sc):De((b=+n)*Aa),u()):b},t.clipExtent=function(n){return arguments.length?(w=n,_=n?Ue(n[0][0],n[0][1],n[1][0],n[1][1]):wt,u()):w},t.scale=function(n){return arguments.length?(h=+n,r()):h},t.translate=function(n){return arguments.length?(g=+n[0],p=+n[1],r()):[g,p]},t.center=function(n){return arguments.length?(v=n[0]%360*Aa,d=n[1]%360*Aa,r()):[v*Ca,d*Ca]},t.rotate=function(n){return arguments.length?(m=n[0]%360*Aa,y=n[1]%360*Aa,x=n.length>2?n[2]%360*Aa:0,r()):[m*Ca,y*Ca,x*Ca]},Zo.rebind(t,f,"precision"),function(){return i=n.apply(this,arguments),t.invert=i.invert&&e,r()}}function er(n){return Qe(n,function(t,e){n.point(t*Aa,e*Aa)})}function rr(n,t){return[n,t]}function ur(n,t){return[n>ba?n-wa:-ba>n?n+wa:n,t]}function ir(n,t,e){return n?t||e?je(ar(n),cr(t,e)):ar(n):t||e?cr(t,e):ur}function or(n){return function(t,e){return t+=n,[t>ba?t-wa:-ba>t?t+wa:t,e]}}function ar(n){var t=or(n);return t.invert=or(-n),t}function cr(n,t){function e(n,t){var e=Math.cos(t),a=Math.cos(n)*e,c=Math.sin(n)*e,s=Math.sin(t),l=s*r+a*u;return[Math.atan2(c*i-l*o,a*r-s*u),G(l*i+c*o)]}var r=Math.cos(n),u=Math.sin(n),i=Math.cos(t),o=Math.sin(t);return e.invert=function(n,t){var e=Math.cos(t),a=Math.cos(n)*e,c=Math.sin(n)*e,s=Math.sin(t),l=s*i-c*o;return[Math.atan2(c*i+s*o,a*r+l*u),G(l*r-a*u)]},e}function sr(n,t){var e=Math.cos(n),r=Math.sin(n);return function(u,i,o,a){var c=o*t;null!=u?(u=lr(e,u),i=lr(e,i),(o>0?i>u:u>i)&&(u+=o*wa)):(u=n+o*wa,i=n-.5*c);for(var s,l=u;o>0?l>i:i>l;l-=c)a.point((s=de([e,-r*Math.cos(l),-r*Math.sin(l)]))[0],s[1])}}function lr(n,t){var e=le(t);e[0]-=n,ve(e);var r=J(-e[1]);return((-e[2]<0?-r:r)+2*Math.PI-ka)%(2*Math.PI)}function fr(n,t,e){var r=Zo.range(n,t-ka,e).concat(t);return function(n){return r.map(function(t){return[n,t]})}}function hr(n,t,e){var r=Zo.range(n,t-ka,e).concat(t);return function(n){return r.map(function(t){return[t,n]})}}function gr(n){return n.source}function pr(n){return n.target}function vr(n,t,e,r){var u=Math.cos(t),i=Math.sin(t),o=Math.cos(r),a=Math.sin(r),c=u*Math.cos(n),s=u*Math.sin(n),l=o*Math.cos(e),f=o*Math.sin(e),h=2*Math.asin(Math.sqrt(tt(r-t)+u*o*tt(e-n))),g=1/Math.sin(h),p=h?function(n){var t=Math.sin(n*=h)*g,e=Math.sin(h-n)*g,r=e*c+t*l,u=e*s+t*f,o=e*i+t*a;return[Math.atan2(u,r)*Ca,Math.atan2(o,Math.sqrt(r*r+u*u))*Ca]}:function(){return[n*Ca,t*Ca]};return p.distance=h,p}function dr(){function n(n,u){var i=Math.sin(u*=Aa),o=Math.cos(u),a=ua((n*=Aa)-t),c=Math.cos(a);Dc+=Math.atan2(Math.sqrt((a=o*Math.sin(a))*a+(a=r*i-e*o*c)*a),e*i+r*o*c),t=n,e=i,r=o}var t,e,r;Pc.point=function(u,i){t=u*Aa,e=Math.sin(i*=Aa),r=Math.cos(i),Pc.point=n},Pc.lineEnd=function(){Pc.point=Pc.lineEnd=v}}function mr(n,t){function e(t,e){var r=Math.cos(t),u=Math.cos(e),i=n(r*u);return[i*u*Math.sin(t),i*Math.sin(e)]}return e.invert=function(n,e){var r=Math.sqrt(n*n+e*e),u=t(r),i=Math.sin(u),o=Math.cos(u);return[Math.atan2(n*i,r*o),Math.asin(r&&e*i/r)]},e}function yr(n,t){function e(n,t){o>0?-Sa+ka>t&&(t=-Sa+ka):t>Sa-ka&&(t=Sa-ka);var e=o/Math.pow(u(t),i);return[e*Math.sin(i*n),o-e*Math.cos(i*n)]}var r=Math.cos(n),u=function(n){return Math.tan(ba/4+n/2)},i=n===t?Math.sin(n):Math.log(r/Math.cos(t))/Math.log(u(t)/u(n)),o=r*Math.pow(u(n),i)/i;return i?(e.invert=function(n,t){var e=o-t,r=B(i)*Math.sqrt(n*n+e*e);return[Math.atan2(n,e)/i,2*Math.atan(Math.pow(o/r,1/i))-Sa]},e):Mr}function xr(n,t){function e(n,t){var e=i-t;return[e*Math.sin(u*n),i-e*Math.cos(u*n)]}var r=Math.cos(n),u=n===t?Math.sin(n):(r-Math.cos(t))/(t-n),i=r/u+n;return ua(u)<ka?rr:(e.invert=function(n,t){var e=i-t;return[Math.atan2(n,e)/u,i-B(u)*Math.sqrt(n*n+e*e)]},e)}function Mr(n,t){return[n,Math.log(Math.tan(ba/4+t/2))]}function _r(n){var t,e=nr(n),r=e.scale,u=e.translate,i=e.clipExtent;return e.scale=function(){var n=r.apply(e,arguments);return n===e?t?e.clipExtent(null):e:n},e.translate=function(){var n=u.apply(e,arguments);return n===e?t?e.clipExtent(null):e:n},e.clipExtent=function(n){var o=i.apply(e,arguments);if(o===e){if(t=null==n){var a=ba*r(),c=u();i([[c[0]-a,c[1]-a],[c[0]+a,c[1]+a]])}}else t&&(o=null);return o},e.clipExtent(null)}function br(n,t){return[Math.log(Math.tan(ba/4+t/2)),-n]}function wr(n){return n[0]}function Sr(n){return n[1]}function kr(n){for(var t=n.length,e=[0,1],r=2,u=2;t>u;u++){for(;r>1&&W(n[e[r-2]],n[e[r-1]],n[u])<=0;)--r;e[r++]=u}return e.slice(0,r)}function Er(n,t){return n[0]-t[0]||n[1]-t[1]}function Ar(n,t,e){return(e[0]-t[0])*(n[1]-t[1])<(e[1]-t[1])*(n[0]-t[0])}function Cr(n,t,e,r){var u=n[0],i=e[0],o=t[0]-u,a=r[0]-i,c=n[1],s=e[1],l=t[1]-c,f=r[1]-s,h=(a*(c-s)-f*(u-i))/(f*o-a*l);return[u+h*o,c+h*l]}function Nr(n){var t=n[0],e=n[n.length-1];return!(t[0]-e[0]||t[1]-e[1])}function zr(){Gr(this),this.edge=this.site=this.circle=null}function Lr(n){var t=Bc.pop()||new zr;return t.site=n,t}function Tr(n){Yr(n),Vc.remove(n),Bc.push(n),Gr(n)}function qr(n){var t=n.circle,e=t.x,r=t.cy,u={x:e,y:r},i=n.P,o=n.N,a=[n];Tr(n);for(var c=i;c.circle&&ua(e-c.circle.x)<ka&&ua(r-c.circle.cy)<ka;)i=c.P,a.unshift(c),Tr(c),c=i;a.unshift(c),Yr(c);for(var s=o;s.circle&&ua(e-s.circle.x)<ka&&ua(r-s.circle.cy)<ka;)o=s.N,a.push(s),Tr(s),s=o;a.push(s),Yr(s);var l,f=a.length;for(l=1;f>l;++l)s=a[l],c=a[l-1],Br(s.edge,c.site,s.site,u);c=a[0],s=a[f-1],s.edge=Xr(c.site,s.site,null,u),Or(c),Or(s)}function Rr(n){for(var t,e,r,u,i=n.x,o=n.y,a=Vc._;a;)if(r=Dr(a,o)-i,r>ka)a=a.L;else{if(u=i-Pr(a,o),!(u>ka)){r>-ka?(t=a.P,e=a):u>-ka?(t=a,e=a.N):t=e=a;break}if(!a.R){t=a;break}a=a.R}var c=Lr(n);if(Vc.insert(t,c),t||e){if(t===e)return Yr(t),e=Lr(t.site),Vc.insert(c,e),c.edge=e.edge=Xr(t.site,c.site),Or(t),Or(e),void 0;if(!e)return c.edge=Xr(t.site,c.site),void 0;Yr(t),Yr(e);var s=t.site,l=s.x,f=s.y,h=n.x-l,g=n.y-f,p=e.site,v=p.x-l,d=p.y-f,m=2*(h*d-g*v),y=h*h+g*g,x=v*v+d*d,M={x:(d*y-g*x)/m+l,y:(h*x-v*y)/m+f};Br(e.edge,s,p,M),c.edge=Xr(s,n,null,M),e.edge=Xr(n,p,null,M),Or(t),Or(e)}}function Dr(n,t){var e=n.site,r=e.x,u=e.y,i=u-t;if(!i)return r;var o=n.P;if(!o)return-1/0;e=o.site;var a=e.x,c=e.y,s=c-t;if(!s)return a;var l=a-r,f=1/i-1/s,h=l/s;return f?(-h+Math.sqrt(h*h-2*f*(l*l/(-2*s)-c+s/2+u-i/2)))/f+r:(r+a)/2}function Pr(n,t){var e=n.N;if(e)return Dr(e,t);var r=n.site;return r.y===t?r.x:1/0}function Ur(n){this.site=n,this.edges=[]}function jr(n){for(var t,e,r,u,i,o,a,c,s,l,f=n[0][0],h=n[1][0],g=n[0][1],p=n[1][1],v=Zc,d=v.length;d--;)if(i=v[d],i&&i.prepare())for(a=i.edges,c=a.length,o=0;c>o;)l=a[o].end(),r=l.x,u=l.y,s=a[++o%c].start(),t=s.x,e=s.y,(ua(r-t)>ka||ua(u-e)>ka)&&(a.splice(o,0,new Wr($r(i.site,l,ua(r-f)<ka&&p-u>ka?{x:f,y:ua(t-f)<ka?e:p}:ua(u-p)<ka&&h-r>ka?{x:ua(e-p)<ka?t:h,y:p}:ua(r-h)<ka&&u-g>ka?{x:h,y:ua(t-h)<ka?e:g}:ua(u-g)<ka&&r-f>ka?{x:ua(e-g)<ka?t:f,y:g}:null),i.site,null)),++c)}function Hr(n,t){return t.angle-n.angle}function Fr(){Gr(this),this.x=this.y=this.arc=this.site=this.cy=null}function Or(n){var t=n.P,e=n.N;if(t&&e){var r=t.site,u=n.site,i=e.site;if(r!==i){var o=u.x,a=u.y,c=r.x-o,s=r.y-a,l=i.x-o,f=i.y-a,h=2*(c*f-s*l);if(!(h>=-Ea)){var g=c*c+s*s,p=l*l+f*f,v=(f*g-s*p)/h,d=(c*p-l*g)/h,f=d+a,m=Wc.pop()||new Fr;m.arc=n,m.site=u,m.x=v+o,m.y=f+Math.sqrt(v*v+d*d),m.cy=f,n.circle=m;for(var y=null,x=$c._;x;)if(m.y<x.y||m.y===x.y&&m.x<=x.x){if(!x.L){y=x.P;break}x=x.L}else{if(!x.R){y=x;break}x=x.R}$c.insert(y,m),y||(Xc=m)}}}}function Yr(n){var t=n.circle;t&&(t.P||(Xc=t.N),$c.remove(t),Wc.push(t),Gr(t),n.circle=null)}function Ir(n){for(var t,e=Ic,r=Pe(n[0][0],n[0][1],n[1][0],n[1][1]),u=e.length;u--;)t=e[u],(!Zr(t,n)||!r(t)||ua(t.a.x-t.b.x)<ka&&ua(t.a.y-t.b.y)<ka)&&(t.a=t.b=null,e.splice(u,1))}function Zr(n,t){var e=n.b;if(e)return!0;var r,u,i=n.a,o=t[0][0],a=t[1][0],c=t[0][1],s=t[1][1],l=n.l,f=n.r,h=l.x,g=l.y,p=f.x,v=f.y,d=(h+p)/2,m=(g+v)/2;if(v===g){if(o>d||d>=a)return;if(h>p){if(i){if(i.y>=s)return}else i={x:d,y:c};e={x:d,y:s}}else{if(i){if(i.y<c)return}else i={x:d,y:s};e={x:d,y:c}}}else if(r=(h-p)/(v-g),u=m-r*d,-1>r||r>1)if(h>p){if(i){if(i.y>=s)return}else i={x:(c-u)/r,y:c};e={x:(s-u)/r,y:s}}else{if(i){if(i.y<c)return}else i={x:(s-u)/r,y:s};e={x:(c-u)/r,y:c}}else if(v>g){if(i){if(i.x>=a)return}else i={x:o,y:r*o+u};e={x:a,y:r*a+u}}else{if(i){if(i.x<o)return}else i={x:a,y:r*a+u};e={x:o,y:r*o+u}}return n.a=i,n.b=e,!0}function Vr(n,t){this.l=n,this.r=t,this.a=this.b=null}function Xr(n,t,e,r){var u=new Vr(n,t);return Ic.push(u),e&&Br(u,n,t,e),r&&Br(u,t,n,r),Zc[n.i].edges.push(new Wr(u,n,t)),Zc[t.i].edges.push(new Wr(u,t,n)),u}function $r(n,t,e){var r=new Vr(n,null);return r.a=t,r.b=e,Ic.push(r),r}function Br(n,t,e,r){n.a||n.b?n.l===e?n.b=r:n.a=r:(n.a=r,n.l=t,n.r=e)}function Wr(n,t,e){var r=n.a,u=n.b;this.edge=n,this.site=t,this.angle=e?Math.atan2(e.y-t.y,e.x-t.x):n.l===t?Math.atan2(u.x-r.x,r.y-u.y):Math.atan2(r.x-u.x,u.y-r.y)}function Jr(){this._=null}function Gr(n){n.U=n.C=n.L=n.R=n.P=n.N=null}function Kr(n,t){var e=t,r=t.R,u=e.U;u?u.L===e?u.L=r:u.R=r:n._=r,r.U=u,e.U=r,e.R=r.L,e.R&&(e.R.U=e),r.L=e}function Qr(n,t){var e=t,r=t.L,u=e.U;u?u.L===e?u.L=r:u.R=r:n._=r,r.U=u,e.U=r,e.L=r.R,e.L&&(e.L.U=e),r.R=e}function nu(n){for(;n.L;)n=n.L;return n}function tu(n,t){var e,r,u,i=n.sort(eu).pop();for(Ic=[],Zc=new Array(n.length),Vc=new Jr,$c=new Jr;;)if(u=Xc,i&&(!u||i.y<u.y||i.y===u.y&&i.x<u.x))(i.x!==e||i.y!==r)&&(Zc[i.i]=new Ur(i),Rr(i),e=i.x,r=i.y),i=n.pop();else{if(!u)break;qr(u.arc)}t&&(Ir(t),jr(t));var o={cells:Zc,edges:Ic};return Vc=$c=Ic=Zc=null,o}function eu(n,t){return t.y-n.y||t.x-n.x}function ru(n,t,e){return(n.x-e.x)*(t.y-n.y)-(n.x-t.x)*(e.y-n.y)}function uu(n){return n.x}function iu(n){return n.y}function ou(){return{leaf:!0,nodes:[],point:null,x:null,y:null}}function au(n,t,e,r,u,i){if(!n(t,e,r,u,i)){var o=.5*(e+u),a=.5*(r+i),c=t.nodes;c[0]&&au(n,c[0],e,r,o,a),c[1]&&au(n,c[1],o,r,u,a),c[2]&&au(n,c[2],e,a,o,i),c[3]&&au(n,c[3],o,a,u,i)}}function cu(n,t){n=Zo.rgb(n),t=Zo.rgb(t);var e=n.r,r=n.g,u=n.b,i=t.r-e,o=t.g-r,a=t.b-u;return function(n){return"#"+dt(Math.round(e+i*n))+dt(Math.round(r+o*n))+dt(Math.round(u+a*n))}}function su(n,t){var e,r={},u={};for(e in n)e in t?r[e]=hu(n[e],t[e]):u[e]=n[e];for(e in t)e in n||(u[e]=t[e]);return function(n){for(e in r)u[e]=r[e](n);return u}}function lu(n,t){return t-=n=+n,function(e){return n+t*e}}function fu(n,t){var e,r,u,i=Gc.lastIndex=Kc.lastIndex=0,o=-1,a=[],c=[];for(n+="",t+="";(e=Gc.exec(n))&&(r=Kc.exec(t));)(u=r.index)>i&&(u=t.substring(i,u),a[o]?a[o]+=u:a[++o]=u),(e=e[0])===(r=r[0])?a[o]?a[o]+=r:a[++o]=r:(a[++o]=null,c.push({i:o,x:lu(e,r)})),i=Kc.lastIndex;return i<t.length&&(u=t.substring(i),a[o]?a[o]+=u:a[++o]=u),a.length<2?c[0]?(t=c[0].x,function(n){return t(n)+""}):function(){return t}:(t=c.length,function(n){for(var e,r=0;t>r;++r)a[(e=c[r]).i]=e.x(n);return a.join("")})}function hu(n,t){for(var e,r=Zo.interpolators.length;--r>=0&&!(e=Zo.interpolators[r](n,t)););return e}function gu(n,t){var e,r=[],u=[],i=n.length,o=t.length,a=Math.min(n.length,t.length);for(e=0;a>e;++e)r.push(hu(n[e],t[e]));for(;i>e;++e)u[e]=n[e];for(;o>e;++e)u[e]=t[e];return function(n){for(e=0;a>e;++e)u[e]=r[e](n);return u}}function pu(n){return function(t){return 0>=t?0:t>=1?1:n(t)}}function vu(n){return function(t){return 1-n(1-t)}}function du(n){return function(t){return.5*(.5>t?n(2*t):2-n(2-2*t))}}function mu(n){return n*n}function yu(n){return n*n*n}function xu(n){if(0>=n)return 0;if(n>=1)return 1;var t=n*n,e=t*n;return 4*(.5>n?e:3*(n-t)+e-.75)}function Mu(n){return function(t){return Math.pow(t,n)}}function _u(n){return 1-Math.cos(n*Sa)}function bu(n){return Math.pow(2,10*(n-1))}function wu(n){return 1-Math.sqrt(1-n*n)}function Su(n,t){var e;return arguments.length<2&&(t=.45),arguments.length?e=t/wa*Math.asin(1/n):(n=1,e=t/4),function(r){return 1+n*Math.pow(2,-10*r)*Math.sin((r-e)*wa/t)}}function ku(n){return n||(n=1.70158),function(t){return t*t*((n+1)*t-n)}}function Eu(n){return 1/2.75>n?7.5625*n*n:2/2.75>n?7.5625*(n-=1.5/2.75)*n+.75:2.5/2.75>n?7.5625*(n-=2.25/2.75)*n+.9375:7.5625*(n-=2.625/2.75)*n+.984375}function Au(n,t){n=Zo.hcl(n),t=Zo.hcl(t);var e=n.h,r=n.c,u=n.l,i=t.h-e,o=t.c-r,a=t.l-u;return isNaN(o)&&(o=0,r=isNaN(r)?t.c:r),isNaN(i)?(i=0,e=isNaN(e)?t.h:e):i>180?i-=360:-180>i&&(i+=360),function(n){return ot(e+i*n,r+o*n,u+a*n)+""}}function Cu(n,t){n=Zo.hsl(n),t=Zo.hsl(t);var e=n.h,r=n.s,u=n.l,i=t.h-e,o=t.s-r,a=t.l-u;return isNaN(o)&&(o=0,r=isNaN(r)?t.s:r),isNaN(i)?(i=0,e=isNaN(e)?t.h:e):i>180?i-=360:-180>i&&(i+=360),function(n){return ut(e+i*n,r+o*n,u+a*n)+""}}function Nu(n,t){n=Zo.lab(n),t=Zo.lab(t);var e=n.l,r=n.a,u=n.b,i=t.l-e,o=t.a-r,a=t.b-u;return function(n){return ct(e+i*n,r+o*n,u+a*n)+""}}function zu(n,t){return t-=n,function(e){return Math.round(n+t*e)}}function Lu(n){var t=[n.a,n.b],e=[n.c,n.d],r=qu(t),u=Tu(t,e),i=qu(Ru(e,t,-u))||0;t[0]*e[1]<e[0]*t[1]&&(t[0]*=-1,t[1]*=-1,r*=-1,u*=-1),this.rotate=(r?Math.atan2(t[1],t[0]):Math.atan2(-e[0],e[1]))*Ca,this.translate=[n.e,n.f],this.scale=[r,i],this.skew=i?Math.atan2(u,i)*Ca:0}function Tu(n,t){return n[0]*t[0]+n[1]*t[1]}function qu(n){var t=Math.sqrt(Tu(n,n));return t&&(n[0]/=t,n[1]/=t),t}function Ru(n,t,e){return n[0]+=e*t[0],n[1]+=e*t[1],n}function Du(n,t){var e,r=[],u=[],i=Zo.transform(n),o=Zo.transform(t),a=i.translate,c=o.translate,s=i.rotate,l=o.rotate,f=i.skew,h=o.skew,g=i.scale,p=o.scale;return a[0]!=c[0]||a[1]!=c[1]?(r.push("translate(",null,",",null,")"),u.push({i:1,x:lu(a[0],c[0])},{i:3,x:lu(a[1],c[1])})):c[0]||c[1]?r.push("translate("+c+")"):r.push(""),s!=l?(s-l>180?l+=360:l-s>180&&(s+=360),u.push({i:r.push(r.pop()+"rotate(",null,")")-2,x:lu(s,l)})):l&&r.push(r.pop()+"rotate("+l+")"),f!=h?u.push({i:r.push(r.pop()+"skewX(",null,")")-2,x:lu(f,h)}):h&&r.push(r.pop()+"skewX("+h+")"),g[0]!=p[0]||g[1]!=p[1]?(e=r.push(r.pop()+"scale(",null,",",null,")"),u.push({i:e-4,x:lu(g[0],p[0])},{i:e-2,x:lu(g[1],p[1])})):(1!=p[0]||1!=p[1])&&r.push(r.pop()+"scale("+p+")"),e=u.length,function(n){for(var t,i=-1;++i<e;)r[(t=u[i]).i]=t.x(n);return r.join("")}}function Pu(n,t){return t=t-(n=+n)?1/(t-n):0,function(e){return(e-n)*t}}function Uu(n,t){return t=t-(n=+n)?1/(t-n):0,function(e){return Math.max(0,Math.min(1,(e-n)*t))}}function ju(n){for(var t=n.source,e=n.target,r=Fu(t,e),u=[t];t!==r;)t=t.parent,u.push(t);for(var i=u.length;e!==r;)u.splice(i,0,e),e=e.parent;return u}function Hu(n){for(var t=[],e=n.parent;null!=e;)t.push(n),n=e,e=e.parent;return t.push(n),t}function Fu(n,t){if(n===t)return n;for(var e=Hu(n),r=Hu(t),u=e.pop(),i=r.pop(),o=null;u===i;)o=u,u=e.pop(),i=r.pop();return o}function Ou(n){n.fixed|=2}function Yu(n){n.fixed&=-7}function Iu(n){n.fixed|=4,n.px=n.x,n.py=n.y}function Zu(n){n.fixed&=-5}function Vu(n,t,e){var r=0,u=0;if(n.charge=0,!n.leaf)for(var i,o=n.nodes,a=o.length,c=-1;++c<a;)i=o[c],null!=i&&(Vu(i,t,e),n.charge+=i.charge,r+=i.charge*i.cx,u+=i.charge*i.cy);if(n.point){n.leaf||(n.point.x+=Math.random()-.5,n.point.y+=Math.random()-.5);var s=t*e[n.point.index];n.charge+=n.pointCharge=s,r+=s*n.point.x,u+=s*n.point.y}n.cx=r/n.charge,n.cy=u/n.charge}function Xu(n,t){return Zo.rebind(n,t,"sort","children","value"),n.nodes=n,n.links=Ku,n}function $u(n,t){for(var e=[n];null!=(n=e.pop());)if(t(n),(u=n.children)&&(r=u.length))for(var r,u;--r>=0;)e.push(u[r])}function Bu(n,t){for(var e=[n],r=[];null!=(n=e.pop());)if(r.push(n),(i=n.children)&&(u=i.length))for(var u,i,o=-1;++o<u;)e.push(i[o]);for(;null!=(n=r.pop());)t(n)}function Wu(n){return n.children}function Ju(n){return n.value}function Gu(n,t){return t.value-n.value}function Ku(n){return Zo.merge(n.map(function(n){return(n.children||[]).map(function(t){return{source:n,target:t}})}))}function Qu(n){return n.x}function ni(n){return n.y}function ti(n,t,e){n.y0=t,n.y=e}function ei(n){return Zo.range(n.length)}function ri(n){for(var t=-1,e=n[0].length,r=[];++t<e;)r[t]=0;return r}function ui(n){for(var t,e=1,r=0,u=n[0][1],i=n.length;i>e;++e)(t=n[e][1])>u&&(r=e,u=t);return r}function ii(n){return n.reduce(oi,0)}function oi(n,t){return n+t[1]}function ai(n,t){return ci(n,Math.ceil(Math.log(t.length)/Math.LN2+1))}function ci(n,t){for(var e=-1,r=+n[0],u=(n[1]-r)/t,i=[];++e<=t;)i[e]=u*e+r;return i}function si(n){return[Zo.min(n),Zo.max(n)]}function li(n,t){return n.value-t.value}function fi(n,t){var e=n._pack_next;n._pack_next=t,t._pack_prev=n,t._pack_next=e,e._pack_prev=t}function hi(n,t){n._pack_next=t,t._pack_prev=n}function gi(n,t){var e=t.x-n.x,r=t.y-n.y,u=n.r+t.r;return.999*u*u>e*e+r*r}function pi(n){function t(n){l=Math.min(n.x-n.r,l),f=Math.max(n.x+n.r,f),h=Math.min(n.y-n.r,h),g=Math.max(n.y+n.r,g)}if((e=n.children)&&(s=e.length)){var e,r,u,i,o,a,c,s,l=1/0,f=-1/0,h=1/0,g=-1/0;if(e.forEach(vi),r=e[0],r.x=-r.r,r.y=0,t(r),s>1&&(u=e[1],u.x=u.r,u.y=0,t(u),s>2))for(i=e[2],yi(r,u,i),t(i),fi(r,i),r._pack_prev=i,fi(i,u),u=r._pack_next,o=3;s>o;o++){yi(r,u,i=e[o]);var p=0,v=1,d=1;for(a=u._pack_next;a!==u;a=a._pack_next,v++)if(gi(a,i)){p=1;break}if(1==p)for(c=r._pack_prev;c!==a._pack_prev&&!gi(c,i);c=c._pack_prev,d++);p?(d>v||v==d&&u.r<r.r?hi(r,u=a):hi(r=c,u),o--):(fi(r,i),u=i,t(i))}var m=(l+f)/2,y=(h+g)/2,x=0;for(o=0;s>o;o++)i=e[o],i.x-=m,i.y-=y,x=Math.max(x,i.r+Math.sqrt(i.x*i.x+i.y*i.y));n.r=x,e.forEach(di)}}function vi(n){n._pack_next=n._pack_prev=n}function di(n){delete n._pack_next,delete n._pack_prev}function mi(n,t,e,r){var u=n.children;if(n.x=t+=r*n.x,n.y=e+=r*n.y,n.r*=r,u)for(var i=-1,o=u.length;++i<o;)mi(u[i],t,e,r)}function yi(n,t,e){var r=n.r+e.r,u=t.x-n.x,i=t.y-n.y;if(r&&(u||i)){var o=t.r+e.r,a=u*u+i*i;o*=o,r*=r;var c=.5+(r-o)/(2*a),s=Math.sqrt(Math.max(0,2*o*(r+a)-(r-=a)*r-o*o))/(2*a);e.x=n.x+c*u+s*i,e.y=n.y+c*i-s*u}else e.x=n.x+r,e.y=n.y}function xi(n,t){return n.parent==t.parent?1:2}function Mi(n){var t=n.children;return t.length?t[0]:n.t}function _i(n){var t,e=n.children;return(t=e.length)?e[t-1]:n.t}function bi(n,t,e){var r=e/(t.i-n.i);t.c-=r,t.s+=e,n.c+=r,t.z+=e,t.m+=e}function wi(n){for(var t,e=0,r=0,u=n.children,i=u.length;--i>=0;)t=u[i],t.z+=e,t.m+=e,e+=t.s+(r+=t.c)}function Si(n,t,e){return n.a.parent===t.parent?n.a:e}function ki(n){return 1+Zo.max(n,function(n){return n.y})}function Ei(n){return n.reduce(function(n,t){return n+t.x},0)/n.length}function Ai(n){var t=n.children;return t&&t.length?Ai(t[0]):n}function Ci(n){var t,e=n.children;return e&&(t=e.length)?Ci(e[t-1]):n}function Ni(n){return{x:n.x,y:n.y,dx:n.dx,dy:n.dy}}function zi(n,t){var e=n.x+t[3],r=n.y+t[0],u=n.dx-t[1]-t[3],i=n.dy-t[0]-t[2];return 0>u&&(e+=u/2,u=0),0>i&&(r+=i/2,i=0),{x:e,y:r,dx:u,dy:i}}function Li(n){var t=n[0],e=n[n.length-1];return e>t?[t,e]:[e,t]}function Ti(n){return n.rangeExtent?n.rangeExtent():Li(n.range())}function qi(n,t,e,r){var u=e(n[0],n[1]),i=r(t[0],t[1]);return function(n){return i(u(n))}}function Ri(n,t){var e,r=0,u=n.length-1,i=n[r],o=n[u];return i>o&&(e=r,r=u,u=e,e=i,i=o,o=e),n[r]=t.floor(i),n[u]=t.ceil(o),n}function Di(n){return n?{floor:function(t){return Math.floor(t/n)*n},ceil:function(t){return Math.ceil(t/n)*n}}:ss}function Pi(n,t,e,r){var u=[],i=[],o=0,a=Math.min(n.length,t.length)-1;for(n[a]<n[0]&&(n=n.slice().reverse(),t=t.slice().reverse());++o<=a;)u.push(e(n[o-1],n[o])),i.push(r(t[o-1],t[o]));return function(t){var e=Zo.bisect(n,t,1,a)-1;return i[e](u[e](t))}}function Ui(n,t,e,r){function u(){var u=Math.min(n.length,t.length)>2?Pi:qi,c=r?Uu:Pu;return o=u(n,t,c,e),a=u(t,n,c,hu),i}function i(n){return o(n)}var o,a;return i.invert=function(n){return a(n)},i.domain=function(t){return arguments.length?(n=t.map(Number),u()):n},i.range=function(n){return arguments.length?(t=n,u()):t},i.rangeRound=function(n){return i.range(n).interpolate(zu)},i.clamp=function(n){return arguments.length?(r=n,u()):r},i.interpolate=function(n){return arguments.length?(e=n,u()):e},i.ticks=function(t){return Oi(n,t)},i.tickFormat=function(t,e){return Yi(n,t,e)},i.nice=function(t){return Hi(n,t),u()},i.copy=function(){return Ui(n,t,e,r)},u()}function ji(n,t){return Zo.rebind(n,t,"range","rangeRound","interpolate","clamp")}function Hi(n,t){return Ri(n,Di(Fi(n,t)[2]))}function Fi(n,t){null==t&&(t=10);var e=Li(n),r=e[1]-e[0],u=Math.pow(10,Math.floor(Math.log(r/t)/Math.LN10)),i=t/r*u;return.15>=i?u*=10:.35>=i?u*=5:.75>=i&&(u*=2),e[0]=Math.ceil(e[0]/u)*u,e[1]=Math.floor(e[1]/u)*u+.5*u,e[2]=u,e}function Oi(n,t){return Zo.range.apply(Zo,Fi(n,t))}function Yi(n,t,e){var r=Fi(n,t);if(e){var u=Ga.exec(e);if(u.shift(),"s"===u[8]){var i=Zo.formatPrefix(Math.max(ua(r[0]),ua(r[1])));return u[7]||(u[7]="."+Ii(i.scale(r[2]))),u[8]="f",e=Zo.format(u.join("")),function(n){return e(i.scale(n))+i.symbol}}u[7]||(u[7]="."+Zi(u[8],r)),e=u.join("")}else e=",."+Ii(r[2])+"f";return Zo.format(e)}function Ii(n){return-Math.floor(Math.log(n)/Math.LN10+.01)}function Zi(n,t){var e=Ii(t[2]);return n in ls?Math.abs(e-Ii(Math.max(ua(t[0]),ua(t[1]))))+ +("e"!==n):e-2*("%"===n)}function Vi(n,t,e,r){function u(n){return(e?Math.log(0>n?0:n):-Math.log(n>0?0:-n))/Math.log(t)}function i(n){return e?Math.pow(t,n):-Math.pow(t,-n)}function o(t){return n(u(t))}return o.invert=function(t){return i(n.invert(t))},o.domain=function(t){return arguments.length?(e=t[0]>=0,n.domain((r=t.map(Number)).map(u)),o):r},o.base=function(e){return arguments.length?(t=+e,n.domain(r.map(u)),o):t},o.nice=function(){var t=Ri(r.map(u),e?Math:hs);return n.domain(t),r=t.map(i),o},o.ticks=function(){var n=Li(r),o=[],a=n[0],c=n[1],s=Math.floor(u(a)),l=Math.ceil(u(c)),f=t%1?2:t;if(isFinite(l-s)){if(e){for(;l>s;s++)for(var h=1;f>h;h++)o.push(i(s)*h);o.push(i(s))}else for(o.push(i(s));s++<l;)for(var h=f-1;h>0;h--)o.push(i(s)*h);for(s=0;o[s]<a;s++);for(l=o.length;o[l-1]>c;l--);o=o.slice(s,l)}return o},o.tickFormat=function(n,t){if(!arguments.length)return fs;arguments.length<2?t=fs:"function"!=typeof t&&(t=Zo.format(t));var r,a=Math.max(.1,n/o.ticks().length),c=e?(r=1e-12,Math.ceil):(r=-1e-12,Math.floor);return function(n){return n/i(c(u(n)+r))<=a?t(n):""}},o.copy=function(){return Vi(n.copy(),t,e,r)},ji(o,n)}function Xi(n,t,e){function r(t){return n(u(t))}var u=$i(t),i=$i(1/t);return r.invert=function(t){return i(n.invert(t))},r.domain=function(t){return arguments.length?(n.domain((e=t.map(Number)).map(u)),r):e},r.ticks=function(n){return Oi(e,n)},r.tickFormat=function(n,t){return Yi(e,n,t)},r.nice=function(n){return r.domain(Hi(e,n))},r.exponent=function(o){return arguments.length?(u=$i(t=o),i=$i(1/t),n.domain(e.map(u)),r):t},r.copy=function(){return Xi(n.copy(),t,e)},ji(r,n)}function $i(n){return function(t){return 0>t?-Math.pow(-t,n):Math.pow(t,n)}}function Bi(n,t){function e(e){return i[((u.get(e)||("range"===t.t?u.set(e,n.push(e)):0/0))-1)%i.length]}function r(t,e){return Zo.range(n.length).map(function(n){return t+e*n})}var u,i,a;return e.domain=function(r){if(!arguments.length)return n;n=[],u=new o;for(var i,a=-1,c=r.length;++a<c;)u.has(i=r[a])||u.set(i,n.push(i));return e[t.t].apply(e,t.a)},e.range=function(n){return arguments.length?(i=n,a=0,t={t:"range",a:arguments},e):i},e.rangePoints=function(u,o){arguments.length<2&&(o=0);var c=u[0],s=u[1],l=(s-c)/(Math.max(1,n.length-1)+o);return i=r(n.length<2?(c+s)/2:c+l*o/2,l),a=0,t={t:"rangePoints",a:arguments},e},e.rangeBands=function(u,o,c){arguments.length<2&&(o=0),arguments.length<3&&(c=o);var s=u[1]<u[0],l=u[s-0],f=u[1-s],h=(f-l)/(n.length-o+2*c);return i=r(l+h*c,h),s&&i.reverse(),a=h*(1-o),t={t:"rangeBands",a:arguments},e},e.rangeRoundBands=function(u,o,c){arguments.length<2&&(o=0),arguments.length<3&&(c=o);var s=u[1]<u[0],l=u[s-0],f=u[1-s],h=Math.floor((f-l)/(n.length-o+2*c)),g=f-l-(n.length-o)*h;return i=r(l+Math.round(g/2),h),s&&i.reverse(),a=Math.round(h*(1-o)),t={t:"rangeRoundBands",a:arguments},e},e.rangeBand=function(){return a},e.rangeExtent=function(){return Li(t.a[0])},e.copy=function(){return Bi(n,t)},e.domain(n)}function Wi(e,r){function u(){var n=0,t=r.length;for(o=[];++n<t;)o[n-1]=Zo.quantile(e,n/t);return i}function i(n){return isNaN(n=+n)?void 0:r[Zo.bisect(o,n)]}var o;return i.domain=function(r){return arguments.length?(e=r.filter(t).sort(n),u()):e},i.range=function(n){return arguments.length?(r=n,u()):r},i.quantiles=function(){return o},i.invertExtent=function(n){return n=r.indexOf(n),0>n?[0/0,0/0]:[n>0?o[n-1]:e[0],n<o.length?o[n]:e[e.length-1]]},i.copy=function(){return Wi(e,r)},u()}function Ji(n,t,e){function r(t){return e[Math.max(0,Math.min(o,Math.floor(i*(t-n))))]}function u(){return i=e.length/(t-n),o=e.length-1,r}var i,o;return r.domain=function(e){return arguments.length?(n=+e[0],t=+e[e.length-1],u()):[n,t]},r.range=function(n){return arguments.length?(e=n,u()):e},r.invertExtent=function(t){return t=e.indexOf(t),t=0>t?0/0:t/i+n,[t,t+1/i]},r.copy=function(){return Ji(n,t,e)},u()}function Gi(n,t){function e(e){return e>=e?t[Zo.bisect(n,e)]:void 0}return e.domain=function(t){return arguments.length?(n=t,e):n},e.range=function(n){return arguments.length?(t=n,e):t},e.invertExtent=function(e){return e=t.indexOf(e),[n[e-1],n[e]]},e.copy=function(){return Gi(n,t)},e}function Ki(n){function t(n){return+n}return t.invert=t,t.domain=t.range=function(e){return arguments.length?(n=e.map(t),t):n},t.ticks=function(t){return Oi(n,t)},t.tickFormat=function(t,e){return Yi(n,t,e)},t.copy=function(){return Ki(n)},t}function Qi(n){return n.innerRadius}function no(n){return n.outerRadius}function to(n){return n.startAngle}function eo(n){return n.endAngle}function ro(n){function t(t){function o(){s.push("M",i(n(l),a))}for(var c,s=[],l=[],f=-1,h=t.length,g=bt(e),p=bt(r);++f<h;)u.call(this,c=t[f],f)?l.push([+g.call(this,c,f),+p.call(this,c,f)]):l.length&&(o(),l=[]);return l.length&&o(),s.length?s.join(""):null}var e=wr,r=Sr,u=we,i=uo,o=i.key,a=.7;return t.x=function(n){return arguments.length?(e=n,t):e},t.y=function(n){return arguments.length?(r=n,t):r},t.defined=function(n){return arguments.length?(u=n,t):u},t.interpolate=function(n){return arguments.length?(o="function"==typeof n?i=n:(i=xs.get(n)||uo).key,t):o},t.tension=function(n){return arguments.length?(a=n,t):a},t}function uo(n){return n.join("L")}function io(n){return uo(n)+"Z"}function oo(n){for(var t=0,e=n.length,r=n[0],u=[r[0],",",r[1]];++t<e;)u.push("H",(r[0]+(r=n[t])[0])/2,"V",r[1]);return e>1&&u.push("H",r[0]),u.join("")}function ao(n){for(var t=0,e=n.length,r=n[0],u=[r[0],",",r[1]];++t<e;)u.push("V",(r=n[t])[1],"H",r[0]);return u.join("")}function co(n){for(var t=0,e=n.length,r=n[0],u=[r[0],",",r[1]];++t<e;)u.push("H",(r=n[t])[0],"V",r[1]);return u.join("")}function so(n,t){return n.length<4?uo(n):n[1]+ho(n.slice(1,n.length-1),go(n,t))}function lo(n,t){return n.length<3?uo(n):n[0]+ho((n.push(n[0]),n),go([n[n.length-2]].concat(n,[n[1]]),t))}function fo(n,t){return n.length<3?uo(n):n[0]+ho(n,go(n,t))}function ho(n,t){if(t.length<1||n.length!=t.length&&n.length!=t.length+2)return uo(n);var e=n.length!=t.length,r="",u=n[0],i=n[1],o=t[0],a=o,c=1;if(e&&(r+="Q"+(i[0]-2*o[0]/3)+","+(i[1]-2*o[1]/3)+","+i[0]+","+i[1],u=n[1],c=2),t.length>1){a=t[1],i=n[c],c++,r+="C"+(u[0]+o[0])+","+(u[1]+o[1])+","+(i[0]-a[0])+","+(i[1]-a[1])+","+i[0]+","+i[1];for(var s=2;s<t.length;s++,c++)i=n[c],a=t[s],r+="S"+(i[0]-a[0])+","+(i[1]-a[1])+","+i[0]+","+i[1]}if(e){var l=n[c];r+="Q"+(i[0]+2*a[0]/3)+","+(i[1]+2*a[1]/3)+","+l[0]+","+l[1]}return r}function go(n,t){for(var e,r=[],u=(1-t)/2,i=n[0],o=n[1],a=1,c=n.length;++a<c;)e=i,i=o,o=n[a],r.push([u*(o[0]-e[0]),u*(o[1]-e[1])]);return r}function po(n){if(n.length<3)return uo(n);var t=1,e=n.length,r=n[0],u=r[0],i=r[1],o=[u,u,u,(r=n[1])[0]],a=[i,i,i,r[1]],c=[u,",",i,"L",xo(bs,o),",",xo(bs,a)];for(n.push(n[e-1]);++t<=e;)r=n[t],o.shift(),o.push(r[0]),a.shift(),a.push(r[1]),Mo(c,o,a);return n.pop(),c.push("L",r),c.join("")}function vo(n){if(n.length<4)return uo(n);for(var t,e=[],r=-1,u=n.length,i=[0],o=[0];++r<3;)t=n[r],i.push(t[0]),o.push(t[1]);for(e.push(xo(bs,i)+","+xo(bs,o)),--r;++r<u;)t=n[r],i.shift(),i.push(t[0]),o.shift(),o.push(t[1]),Mo(e,i,o);return e.join("")}function mo(n){for(var t,e,r=-1,u=n.length,i=u+4,o=[],a=[];++r<4;)e=n[r%u],o.push(e[0]),a.push(e[1]);for(t=[xo(bs,o),",",xo(bs,a)],--r;++r<i;)e=n[r%u],o.shift(),o.push(e[0]),a.shift(),a.push(e[1]),Mo(t,o,a);return t.join("")}function yo(n,t){var e=n.length-1;if(e)for(var r,u,i=n[0][0],o=n[0][1],a=n[e][0]-i,c=n[e][1]-o,s=-1;++s<=e;)r=n[s],u=s/e,r[0]=t*r[0]+(1-t)*(i+u*a),r[1]=t*r[1]+(1-t)*(o+u*c);return po(n)}function xo(n,t){return n[0]*t[0]+n[1]*t[1]+n[2]*t[2]+n[3]*t[3]}function Mo(n,t,e){n.push("C",xo(Ms,t),",",xo(Ms,e),",",xo(_s,t),",",xo(_s,e),",",xo(bs,t),",",xo(bs,e))}function _o(n,t){return(t[1]-n[1])/(t[0]-n[0])}function bo(n){for(var t=0,e=n.length-1,r=[],u=n[0],i=n[1],o=r[0]=_o(u,i);++t<e;)r[t]=(o+(o=_o(u=i,i=n[t+1])))/2;return r[t]=o,r}function wo(n){for(var t,e,r,u,i=[],o=bo(n),a=-1,c=n.length-1;++a<c;)t=_o(n[a],n[a+1]),ua(t)<ka?o[a]=o[a+1]=0:(e=o[a]/t,r=o[a+1]/t,u=e*e+r*r,u>9&&(u=3*t/Math.sqrt(u),o[a]=u*e,o[a+1]=u*r));for(a=-1;++a<=c;)u=(n[Math.min(c,a+1)][0]-n[Math.max(0,a-1)][0])/(6*(1+o[a]*o[a])),i.push([u||0,o[a]*u||0]);return i}function So(n){return n.length<3?uo(n):n[0]+ho(n,wo(n))}function ko(n){for(var t,e,r,u=-1,i=n.length;++u<i;)t=n[u],e=t[0],r=t[1]+ms,t[0]=e*Math.cos(r),t[1]=e*Math.sin(r);return n}function Eo(n){function t(t){function c(){v.push("M",a(n(m),f),l,s(n(d.reverse()),f),"Z")}for(var h,g,p,v=[],d=[],m=[],y=-1,x=t.length,M=bt(e),_=bt(u),b=e===r?function(){return g}:bt(r),w=u===i?function(){return p}:bt(i);++y<x;)o.call(this,h=t[y],y)?(d.push([g=+M.call(this,h,y),p=+_.call(this,h,y)]),m.push([+b.call(this,h,y),+w.call(this,h,y)])):d.length&&(c(),d=[],m=[]);return d.length&&c(),v.length?v.join(""):null}var e=wr,r=wr,u=0,i=Sr,o=we,a=uo,c=a.key,s=a,l="L",f=.7;return t.x=function(n){return arguments.length?(e=r=n,t):r},t.x0=function(n){return arguments.length?(e=n,t):e},t.x1=function(n){return arguments.length?(r=n,t):r},t.y=function(n){return arguments.length?(u=i=n,t):i},t.y0=function(n){return arguments.length?(u=n,t):u},t.y1=function(n){return arguments.length?(i=n,t):i},t.defined=function(n){return arguments.length?(o=n,t):o},t.interpolate=function(n){return arguments.length?(c="function"==typeof n?a=n:(a=xs.get(n)||uo).key,s=a.reverse||a,l=a.closed?"M":"L",t):c},t.tension=function(n){return arguments.length?(f=n,t):f},t}function Ao(n){return n.radius}function Co(n){return[n.x,n.y]}function No(n){return function(){var t=n.apply(this,arguments),e=t[0],r=t[1]+ms;return[e*Math.cos(r),e*Math.sin(r)]}}function zo(){return 64}function Lo(){return"circle"}function To(n){var t=Math.sqrt(n/ba);return"M0,"+t+"A"+t+","+t+" 0 1,1 0,"+-t+"A"+t+","+t+" 0 1,1 0,"+t+"Z"}function qo(n,t){return sa(n,Cs),n.id=t,n}function Ro(n,t,e,r){var u=n.id;return P(n,"function"==typeof e?function(n,i,o){n.__transition__[u].tween.set(t,r(e.call(n,n.__data__,i,o)))}:(e=r(e),function(n){n.__transition__[u].tween.set(t,e)}))}function Do(n){return null==n&&(n=""),function(){this.textContent=n}}function Po(n,t,e,r){var u=n.__transition__||(n.__transition__={active:0,count:0}),i=u[e];if(!i){var a=r.time;i=u[e]={tween:new o,time:a,ease:r.ease,delay:r.delay,duration:r.duration},++u.count,Zo.timer(function(r){function o(r){return u.active>e?s():(u.active=e,i.event&&i.event.start.call(n,l,t),i.tween.forEach(function(e,r){(r=r.call(n,l,t))&&v.push(r)}),Zo.timer(function(){return p.c=c(r||1)?we:c,1},0,a),void 0)}function c(r){if(u.active!==e)return s();for(var o=r/g,a=f(o),c=v.length;c>0;)v[--c].call(n,a);
return o>=1?(i.event&&i.event.end.call(n,l,t),s()):void 0}function s(){return--u.count?delete u[e]:delete n.__transition__,1}var l=n.__data__,f=i.ease,h=i.delay,g=i.duration,p=Ba,v=[];return p.t=h+a,r>=h?o(r-h):(p.c=o,void 0)},0,a)}}function Uo(n,t){n.attr("transform",function(n){return"translate("+t(n)+",0)"})}function jo(n,t){n.attr("transform",function(n){return"translate(0,"+t(n)+")"})}function Ho(n){return n.toISOString()}function Fo(n,t,e){function r(t){return n(t)}function u(n,e){var r=n[1]-n[0],u=r/e,i=Zo.bisect(Us,u);return i==Us.length?[t.year,Fi(n.map(function(n){return n/31536e6}),e)[2]]:i?t[u/Us[i-1]<Us[i]/u?i-1:i]:[Fs,Fi(n,e)[2]]}return r.invert=function(t){return Oo(n.invert(t))},r.domain=function(t){return arguments.length?(n.domain(t),r):n.domain().map(Oo)},r.nice=function(n,t){function e(e){return!isNaN(e)&&!n.range(e,Oo(+e+1),t).length}var i=r.domain(),o=Li(i),a=null==n?u(o,10):"number"==typeof n&&u(o,n);return a&&(n=a[0],t=a[1]),r.domain(Ri(i,t>1?{floor:function(t){for(;e(t=n.floor(t));)t=Oo(t-1);return t},ceil:function(t){for(;e(t=n.ceil(t));)t=Oo(+t+1);return t}}:n))},r.ticks=function(n,t){var e=Li(r.domain()),i=null==n?u(e,10):"number"==typeof n?u(e,n):!n.range&&[{range:n},t];return i&&(n=i[0],t=i[1]),n.range(e[0],Oo(+e[1]+1),1>t?1:t)},r.tickFormat=function(){return e},r.copy=function(){return Fo(n.copy(),t,e)},ji(r,n)}function Oo(n){return new Date(n)}function Yo(n){return JSON.parse(n.responseText)}function Io(n){var t=$o.createRange();return t.selectNode($o.body),t.createContextualFragment(n.responseText)}var Zo={version:"3.4.11"};Date.now||(Date.now=function(){return+new Date});var Vo=[].slice,Xo=function(n){return Vo.call(n)},$o=document,Bo=$o.documentElement,Wo=window;try{Xo(Bo.childNodes)[0].nodeType}catch(Jo){Xo=function(n){for(var t=n.length,e=new Array(t);t--;)e[t]=n[t];return e}}try{$o.createElement("div").style.setProperty("opacity",0,"")}catch(Go){var Ko=Wo.Element.prototype,Qo=Ko.setAttribute,na=Ko.setAttributeNS,ta=Wo.CSSStyleDeclaration.prototype,ea=ta.setProperty;Ko.setAttribute=function(n,t){Qo.call(this,n,t+"")},Ko.setAttributeNS=function(n,t,e){na.call(this,n,t,e+"")},ta.setProperty=function(n,t,e){ea.call(this,n,t+"",e)}}Zo.ascending=n,Zo.descending=function(n,t){return n>t?-1:t>n?1:t>=n?0:0/0},Zo.min=function(n,t){var e,r,u=-1,i=n.length;if(1===arguments.length){for(;++u<i&&!(null!=(e=n[u])&&e>=e);)e=void 0;for(;++u<i;)null!=(r=n[u])&&e>r&&(e=r)}else{for(;++u<i&&!(null!=(e=t.call(n,n[u],u))&&e>=e);)e=void 0;for(;++u<i;)null!=(r=t.call(n,n[u],u))&&e>r&&(e=r)}return e},Zo.max=function(n,t){var e,r,u=-1,i=n.length;if(1===arguments.length){for(;++u<i&&!(null!=(e=n[u])&&e>=e);)e=void 0;for(;++u<i;)null!=(r=n[u])&&r>e&&(e=r)}else{for(;++u<i&&!(null!=(e=t.call(n,n[u],u))&&e>=e);)e=void 0;for(;++u<i;)null!=(r=t.call(n,n[u],u))&&r>e&&(e=r)}return e},Zo.extent=function(n,t){var e,r,u,i=-1,o=n.length;if(1===arguments.length){for(;++i<o&&!(null!=(e=u=n[i])&&e>=e);)e=u=void 0;for(;++i<o;)null!=(r=n[i])&&(e>r&&(e=r),r>u&&(u=r))}else{for(;++i<o&&!(null!=(e=u=t.call(n,n[i],i))&&e>=e);)e=void 0;for(;++i<o;)null!=(r=t.call(n,n[i],i))&&(e>r&&(e=r),r>u&&(u=r))}return[e,u]},Zo.sum=function(n,t){var e,r=0,u=n.length,i=-1;if(1===arguments.length)for(;++i<u;)isNaN(e=+n[i])||(r+=e);else for(;++i<u;)isNaN(e=+t.call(n,n[i],i))||(r+=e);return r},Zo.mean=function(n,e){var r,u=0,i=n.length,o=-1,a=i;if(1===arguments.length)for(;++o<i;)t(r=n[o])?u+=r:--a;else for(;++o<i;)t(r=e.call(n,n[o],o))?u+=r:--a;return a?u/a:void 0},Zo.quantile=function(n,t){var e=(n.length-1)*t+1,r=Math.floor(e),u=+n[r-1],i=e-r;return i?u+i*(n[r]-u):u},Zo.median=function(e,r){return arguments.length>1&&(e=e.map(r)),e=e.filter(t),e.length?Zo.quantile(e.sort(n),.5):void 0};var ra=e(n);Zo.bisectLeft=ra.left,Zo.bisect=Zo.bisectRight=ra.right,Zo.bisector=function(t){return e(1===t.length?function(e,r){return n(t(e),r)}:t)},Zo.shuffle=function(n){for(var t,e,r=n.length;r;)e=0|Math.random()*r--,t=n[r],n[r]=n[e],n[e]=t;return n},Zo.permute=function(n,t){for(var e=t.length,r=new Array(e);e--;)r[e]=n[t[e]];return r},Zo.pairs=function(n){for(var t,e=0,r=n.length-1,u=n[0],i=new Array(0>r?0:r);r>e;)i[e]=[t=u,u=n[++e]];return i},Zo.zip=function(){if(!(u=arguments.length))return[];for(var n=-1,t=Zo.min(arguments,r),e=new Array(t);++n<t;)for(var u,i=-1,o=e[n]=new Array(u);++i<u;)o[i]=arguments[i][n];return e},Zo.transpose=function(n){return Zo.zip.apply(Zo,n)},Zo.keys=function(n){var t=[];for(var e in n)t.push(e);return t},Zo.values=function(n){var t=[];for(var e in n)t.push(n[e]);return t},Zo.entries=function(n){var t=[];for(var e in n)t.push({key:e,value:n[e]});return t},Zo.merge=function(n){for(var t,e,r,u=n.length,i=-1,o=0;++i<u;)o+=n[i].length;for(e=new Array(o);--u>=0;)for(r=n[u],t=r.length;--t>=0;)e[--o]=r[t];return e};var ua=Math.abs;Zo.range=function(n,t,e){if(arguments.length<3&&(e=1,arguments.length<2&&(t=n,n=0)),1/0===(t-n)/e)throw new Error("infinite range");var r,i=[],o=u(ua(e)),a=-1;if(n*=o,t*=o,e*=o,0>e)for(;(r=n+e*++a)>t;)i.push(r/o);else for(;(r=n+e*++a)<t;)i.push(r/o);return i},Zo.map=function(n){var t=new o;if(n instanceof o)n.forEach(function(n,e){t.set(n,e)});else for(var e in n)t.set(e,n[e]);return t},i(o,{has:a,get:function(n){return this[ia+n]},set:function(n,t){return this[ia+n]=t},remove:c,keys:s,values:function(){var n=[];return this.forEach(function(t,e){n.push(e)}),n},entries:function(){var n=[];return this.forEach(function(t,e){n.push({key:t,value:e})}),n},size:l,empty:f,forEach:function(n){for(var t in this)t.charCodeAt(0)===oa&&n.call(this,t.substring(1),this[t])}});var ia="\x00",oa=ia.charCodeAt(0);Zo.nest=function(){function n(t,a,c){if(c>=i.length)return r?r.call(u,a):e?a.sort(e):a;for(var s,l,f,h,g=-1,p=a.length,v=i[c++],d=new o;++g<p;)(h=d.get(s=v(l=a[g])))?h.push(l):d.set(s,[l]);return t?(l=t(),f=function(e,r){l.set(e,n(t,r,c))}):(l={},f=function(e,r){l[e]=n(t,r,c)}),d.forEach(f),l}function t(n,e){if(e>=i.length)return n;var r=[],u=a[e++];return n.forEach(function(n,u){r.push({key:n,values:t(u,e)})}),u?r.sort(function(n,t){return u(n.key,t.key)}):r}var e,r,u={},i=[],a=[];return u.map=function(t,e){return n(e,t,0)},u.entries=function(e){return t(n(Zo.map,e,0),0)},u.key=function(n){return i.push(n),u},u.sortKeys=function(n){return a[i.length-1]=n,u},u.sortValues=function(n){return e=n,u},u.rollup=function(n){return r=n,u},u},Zo.set=function(n){var t=new h;if(n)for(var e=0,r=n.length;r>e;++e)t.add(n[e]);return t},i(h,{has:a,add:function(n){return this[ia+n]=!0,n},remove:function(n){return n=ia+n,n in this&&delete this[n]},values:s,size:l,empty:f,forEach:function(n){for(var t in this)t.charCodeAt(0)===oa&&n.call(this,t.substring(1))}}),Zo.behavior={},Zo.rebind=function(n,t){for(var e,r=1,u=arguments.length;++r<u;)n[e=arguments[r]]=g(n,t,t[e]);return n};var aa=["webkit","ms","moz","Moz","o","O"];Zo.dispatch=function(){for(var n=new d,t=-1,e=arguments.length;++t<e;)n[arguments[t]]=m(n);return n},d.prototype.on=function(n,t){var e=n.indexOf("."),r="";if(e>=0&&(r=n.substring(e+1),n=n.substring(0,e)),n)return arguments.length<2?this[n].on(r):this[n].on(r,t);if(2===arguments.length){if(null==t)for(n in this)this.hasOwnProperty(n)&&this[n].on(r,null);return this}},Zo.event=null,Zo.requote=function(n){return n.replace(ca,"\\$&")};var ca=/[\\\^\$\*\+\?\|\[\]\(\)\.\{\}]/g,sa={}.__proto__?function(n,t){n.__proto__=t}:function(n,t){for(var e in t)n[e]=t[e]},la=function(n,t){return t.querySelector(n)},fa=function(n,t){return t.querySelectorAll(n)},ha=Bo.matches||Bo[p(Bo,"matchesSelector")],ga=function(n,t){return ha.call(n,t)};"function"==typeof Sizzle&&(la=function(n,t){return Sizzle(n,t)[0]||null},fa=Sizzle,ga=Sizzle.matchesSelector),Zo.selection=function(){return ma};var pa=Zo.selection.prototype=[];pa.select=function(n){var t,e,r,u,i=[];n=b(n);for(var o=-1,a=this.length;++o<a;){i.push(t=[]),t.parentNode=(r=this[o]).parentNode;for(var c=-1,s=r.length;++c<s;)(u=r[c])?(t.push(e=n.call(u,u.__data__,c,o)),e&&"__data__"in u&&(e.__data__=u.__data__)):t.push(null)}return _(i)},pa.selectAll=function(n){var t,e,r=[];n=w(n);for(var u=-1,i=this.length;++u<i;)for(var o=this[u],a=-1,c=o.length;++a<c;)(e=o[a])&&(r.push(t=Xo(n.call(e,e.__data__,a,u))),t.parentNode=e);return _(r)};var va={svg:"http://www.w3.org/2000/svg",xhtml:"http://www.w3.org/1999/xhtml",xlink:"http://www.w3.org/1999/xlink",xml:"http://www.w3.org/XML/1998/namespace",xmlns:"http://www.w3.org/2000/xmlns/"};Zo.ns={prefix:va,qualify:function(n){var t=n.indexOf(":"),e=n;return t>=0&&(e=n.substring(0,t),n=n.substring(t+1)),va.hasOwnProperty(e)?{space:va[e],local:n}:n}},pa.attr=function(n,t){if(arguments.length<2){if("string"==typeof n){var e=this.node();return n=Zo.ns.qualify(n),n.local?e.getAttributeNS(n.space,n.local):e.getAttribute(n)}for(t in n)this.each(S(t,n[t]));return this}return this.each(S(n,t))},pa.classed=function(n,t){if(arguments.length<2){if("string"==typeof n){var e=this.node(),r=(n=A(n)).length,u=-1;if(t=e.classList){for(;++u<r;)if(!t.contains(n[u]))return!1}else for(t=e.getAttribute("class");++u<r;)if(!E(n[u]).test(t))return!1;return!0}for(t in n)this.each(C(t,n[t]));return this}return this.each(C(n,t))},pa.style=function(n,t,e){var r=arguments.length;if(3>r){if("string"!=typeof n){2>r&&(t="");for(e in n)this.each(z(e,n[e],t));return this}if(2>r)return Wo.getComputedStyle(this.node(),null).getPropertyValue(n);e=""}return this.each(z(n,t,e))},pa.property=function(n,t){if(arguments.length<2){if("string"==typeof n)return this.node()[n];for(t in n)this.each(L(t,n[t]));return this}return this.each(L(n,t))},pa.text=function(n){return arguments.length?this.each("function"==typeof n?function(){var t=n.apply(this,arguments);this.textContent=null==t?"":t}:null==n?function(){this.textContent=""}:function(){this.textContent=n}):this.node().textContent},pa.html=function(n){return arguments.length?this.each("function"==typeof n?function(){var t=n.apply(this,arguments);this.innerHTML=null==t?"":t}:null==n?function(){this.innerHTML=""}:function(){this.innerHTML=n}):this.node().innerHTML},pa.append=function(n){return n=T(n),this.select(function(){return this.appendChild(n.apply(this,arguments))})},pa.insert=function(n,t){return n=T(n),t=b(t),this.select(function(){return this.insertBefore(n.apply(this,arguments),t.apply(this,arguments)||null)})},pa.remove=function(){return this.each(function(){var n=this.parentNode;n&&n.removeChild(this)})},pa.data=function(n,t){function e(n,e){var r,u,i,a=n.length,f=e.length,h=Math.min(a,f),g=new Array(f),p=new Array(f),v=new Array(a);if(t){var d,m=new o,y=new o,x=[];for(r=-1;++r<a;)d=t.call(u=n[r],u.__data__,r),m.has(d)?v[r]=u:m.set(d,u),x.push(d);for(r=-1;++r<f;)d=t.call(e,i=e[r],r),(u=m.get(d))?(g[r]=u,u.__data__=i):y.has(d)||(p[r]=q(i)),y.set(d,i),m.remove(d);for(r=-1;++r<a;)m.has(x[r])&&(v[r]=n[r])}else{for(r=-1;++r<h;)u=n[r],i=e[r],u?(u.__data__=i,g[r]=u):p[r]=q(i);for(;f>r;++r)p[r]=q(e[r]);for(;a>r;++r)v[r]=n[r]}p.update=g,p.parentNode=g.parentNode=v.parentNode=n.parentNode,c.push(p),s.push(g),l.push(v)}var r,u,i=-1,a=this.length;if(!arguments.length){for(n=new Array(a=(r=this[0]).length);++i<a;)(u=r[i])&&(n[i]=u.__data__);return n}var c=U([]),s=_([]),l=_([]);if("function"==typeof n)for(;++i<a;)e(r=this[i],n.call(r,r.parentNode.__data__,i));else for(;++i<a;)e(r=this[i],n);return s.enter=function(){return c},s.exit=function(){return l},s},pa.datum=function(n){return arguments.length?this.property("__data__",n):this.property("__data__")},pa.filter=function(n){var t,e,r,u=[];"function"!=typeof n&&(n=R(n));for(var i=0,o=this.length;o>i;i++){u.push(t=[]),t.parentNode=(e=this[i]).parentNode;for(var a=0,c=e.length;c>a;a++)(r=e[a])&&n.call(r,r.__data__,a,i)&&t.push(r)}return _(u)},pa.order=function(){for(var n=-1,t=this.length;++n<t;)for(var e,r=this[n],u=r.length-1,i=r[u];--u>=0;)(e=r[u])&&(i&&i!==e.nextSibling&&i.parentNode.insertBefore(e,i),i=e);return this},pa.sort=function(n){n=D.apply(this,arguments);for(var t=-1,e=this.length;++t<e;)this[t].sort(n);return this.order()},pa.each=function(n){return P(this,function(t,e,r){n.call(t,t.__data__,e,r)})},pa.call=function(n){var t=Xo(arguments);return n.apply(t[0]=this,t),this},pa.empty=function(){return!this.node()},pa.node=function(){for(var n=0,t=this.length;t>n;n++)for(var e=this[n],r=0,u=e.length;u>r;r++){var i=e[r];if(i)return i}return null},pa.size=function(){var n=0;return this.each(function(){++n}),n};var da=[];Zo.selection.enter=U,Zo.selection.enter.prototype=da,da.append=pa.append,da.empty=pa.empty,da.node=pa.node,da.call=pa.call,da.size=pa.size,da.select=function(n){for(var t,e,r,u,i,o=[],a=-1,c=this.length;++a<c;){r=(u=this[a]).update,o.push(t=[]),t.parentNode=u.parentNode;for(var s=-1,l=u.length;++s<l;)(i=u[s])?(t.push(r[s]=e=n.call(u.parentNode,i.__data__,s,a)),e.__data__=i.__data__):t.push(null)}return _(o)},da.insert=function(n,t){return arguments.length<2&&(t=j(this)),pa.insert.call(this,n,t)},pa.transition=function(){for(var n,t,e=Ss||++Ns,r=[],u=ks||{time:Date.now(),ease:xu,delay:0,duration:250},i=-1,o=this.length;++i<o;){r.push(n=[]);for(var a=this[i],c=-1,s=a.length;++c<s;)(t=a[c])&&Po(t,c,e,u),n.push(t)}return qo(r,e)},pa.interrupt=function(){return this.each(H)},Zo.select=function(n){var t=["string"==typeof n?la(n,$o):n];return t.parentNode=Bo,_([t])},Zo.selectAll=function(n){var t=Xo("string"==typeof n?fa(n,$o):n);return t.parentNode=Bo,_([t])};var ma=Zo.select(Bo);pa.on=function(n,t,e){var r=arguments.length;if(3>r){if("string"!=typeof n){2>r&&(t=!1);for(e in n)this.each(F(e,n[e],t));return this}if(2>r)return(r=this.node()["__on"+n])&&r._;e=!1}return this.each(F(n,t,e))};var ya=Zo.map({mouseenter:"mouseover",mouseleave:"mouseout"});ya.forEach(function(n){"on"+n in $o&&ya.remove(n)});var xa="onselectstart"in $o?null:p(Bo.style,"userSelect"),Ma=0;Zo.mouse=function(n){return Z(n,x())};var _a=/WebKit/.test(Wo.navigator.userAgent)?-1:0;Zo.touches=function(n,t){return arguments.length<2&&(t=x().touches),t?Xo(t).map(function(t){var e=Z(n,t);return e.identifier=t.identifier,e}):[]},Zo.behavior.drag=function(){function n(){this.on("mousedown.drag",u).on("touchstart.drag",i)}function t(n,t,u,i,o){return function(){function a(){var n,e,r=t(h,v);r&&(n=r[0]-x[0],e=r[1]-x[1],p|=n|e,x=r,g({type:"drag",x:r[0]+s[0],y:r[1]+s[1],dx:n,dy:e}))}function c(){t(h,v)&&(m.on(i+d,null).on(o+d,null),y(p&&Zo.event.target===f),g({type:"dragend"}))}var s,l=this,f=Zo.event.target,h=l.parentNode,g=e.of(l,arguments),p=0,v=n(),d=".drag"+(null==v?"":"-"+v),m=Zo.select(u()).on(i+d,a).on(o+d,c),y=I(),x=t(h,v);r?(s=r.apply(l,arguments),s=[s.x-x[0],s.y-x[1]]):s=[0,0],g({type:"dragstart"})}}var e=M(n,"drag","dragstart","dragend"),r=null,u=t(v,Zo.mouse,$,"mousemove","mouseup"),i=t(V,Zo.touch,X,"touchmove","touchend");return n.origin=function(t){return arguments.length?(r=t,n):r},Zo.rebind(n,e,"on")};var ba=Math.PI,wa=2*ba,Sa=ba/2,ka=1e-6,Ea=ka*ka,Aa=ba/180,Ca=180/ba,Na=Math.SQRT2,za=2,La=4;Zo.interpolateZoom=function(n,t){function e(n){var t=n*y;if(m){var e=Q(v),o=i/(za*h)*(e*nt(Na*t+v)-K(v));return[r+o*s,u+o*l,i*e/Q(Na*t+v)]}return[r+n*s,u+n*l,i*Math.exp(Na*t)]}var r=n[0],u=n[1],i=n[2],o=t[0],a=t[1],c=t[2],s=o-r,l=a-u,f=s*s+l*l,h=Math.sqrt(f),g=(c*c-i*i+La*f)/(2*i*za*h),p=(c*c-i*i-La*f)/(2*c*za*h),v=Math.log(Math.sqrt(g*g+1)-g),d=Math.log(Math.sqrt(p*p+1)-p),m=d-v,y=(m||Math.log(c/i))/Na;return e.duration=1e3*y,e},Zo.behavior.zoom=function(){function n(n){n.on(A,s).on(Ra+".zoom",f).on("dblclick.zoom",h).on(z,l)}function t(n){return[(n[0]-S.x)/S.k,(n[1]-S.y)/S.k]}function e(n){return[n[0]*S.k+S.x,n[1]*S.k+S.y]}function r(n){S.k=Math.max(E[0],Math.min(E[1],n))}function u(n,t){t=e(t),S.x+=n[0]-t[0],S.y+=n[1]-t[1]}function i(){_&&_.domain(x.range().map(function(n){return(n-S.x)/S.k}).map(x.invert)),w&&w.domain(b.range().map(function(n){return(n-S.y)/S.k}).map(b.invert))}function o(n){n({type:"zoomstart"})}function a(n){i(),n({type:"zoom",scale:S.k,translate:[S.x,S.y]})}function c(n){n({type:"zoomend"})}function s(){function n(){l=1,u(Zo.mouse(r),h),a(s)}function e(){f.on(C,null).on(N,null),g(l&&Zo.event.target===i),c(s)}var r=this,i=Zo.event.target,s=L.of(r,arguments),l=0,f=Zo.select(Wo).on(C,n).on(N,e),h=t(Zo.mouse(r)),g=I();H.call(r),o(s)}function l(){function n(){var n=Zo.touches(g);return h=S.k,n.forEach(function(n){n.identifier in v&&(v[n.identifier]=t(n))}),n}function e(){var t=Zo.event.target;Zo.select(t).on(M,i).on(_,f),b.push(t);for(var e=Zo.event.changedTouches,o=0,c=e.length;c>o;++o)v[e[o].identifier]=null;var s=n(),l=Date.now();if(1===s.length){if(500>l-m){var h=s[0],g=v[h.identifier];r(2*S.k),u(h,g),y(),a(p)}m=l}else if(s.length>1){var h=s[0],x=s[1],w=h[0]-x[0],k=h[1]-x[1];d=w*w+k*k}}function i(){for(var n,t,e,i,o=Zo.touches(g),c=0,s=o.length;s>c;++c,i=null)if(e=o[c],i=v[e.identifier]){if(t)break;n=e,t=i}if(i){var l=(l=e[0]-n[0])*l+(l=e[1]-n[1])*l,f=d&&Math.sqrt(l/d);n=[(n[0]+e[0])/2,(n[1]+e[1])/2],t=[(t[0]+i[0])/2,(t[1]+i[1])/2],r(f*h)}m=null,u(n,t),a(p)}function f(){if(Zo.event.touches.length){for(var t=Zo.event.changedTouches,e=0,r=t.length;r>e;++e)delete v[t[e].identifier];for(var u in v)return void n()}Zo.selectAll(b).on(x,null),w.on(A,s).on(z,l),k(),c(p)}var h,g=this,p=L.of(g,arguments),v={},d=0,x=".zoom-"+Zo.event.changedTouches[0].identifier,M="touchmove"+x,_="touchend"+x,b=[],w=Zo.select(g).on(A,null).on(z,e),k=I();H.call(g),e(),o(p)}function f(){var n=L.of(this,arguments);d?clearTimeout(d):(g=t(p=v||Zo.mouse(this)),H.call(this),o(n)),d=setTimeout(function(){d=null,c(n)},50),y(),r(Math.pow(2,.002*Ta())*S.k),u(p,g),a(n)}function h(){var n=L.of(this,arguments),e=Zo.mouse(this),i=t(e),s=Math.log(S.k)/Math.LN2;o(n),r(Math.pow(2,Zo.event.shiftKey?Math.ceil(s)-1:Math.floor(s)+1)),u(e,i),a(n),c(n)}var g,p,v,d,m,x,_,b,w,S={x:0,y:0,k:1},k=[960,500],E=qa,A="mousedown.zoom",C="mousemove.zoom",N="mouseup.zoom",z="touchstart.zoom",L=M(n,"zoomstart","zoom","zoomend");return n.event=function(n){n.each(function(){var n=L.of(this,arguments),t=S;Ss?Zo.select(this).transition().each("start.zoom",function(){S=this.__chart__||{x:0,y:0,k:1},o(n)}).tween("zoom:zoom",function(){var e=k[0],r=k[1],u=e/2,i=r/2,o=Zo.interpolateZoom([(u-S.x)/S.k,(i-S.y)/S.k,e/S.k],[(u-t.x)/t.k,(i-t.y)/t.k,e/t.k]);return function(t){var r=o(t),c=e/r[2];this.__chart__=S={x:u-r[0]*c,y:i-r[1]*c,k:c},a(n)}}).each("end.zoom",function(){c(n)}):(this.__chart__=S,o(n),a(n),c(n))})},n.translate=function(t){return arguments.length?(S={x:+t[0],y:+t[1],k:S.k},i(),n):[S.x,S.y]},n.scale=function(t){return arguments.length?(S={x:S.x,y:S.y,k:+t},i(),n):S.k},n.scaleExtent=function(t){return arguments.length?(E=null==t?qa:[+t[0],+t[1]],n):E},n.center=function(t){return arguments.length?(v=t&&[+t[0],+t[1]],n):v},n.size=function(t){return arguments.length?(k=t&&[+t[0],+t[1]],n):k},n.x=function(t){return arguments.length?(_=t,x=t.copy(),S={x:0,y:0,k:1},n):_},n.y=function(t){return arguments.length?(w=t,b=t.copy(),S={x:0,y:0,k:1},n):w},Zo.rebind(n,L,"on")};var Ta,qa=[0,1/0],Ra="onwheel"in $o?(Ta=function(){return-Zo.event.deltaY*(Zo.event.deltaMode?120:1)},"wheel"):"onmousewheel"in $o?(Ta=function(){return Zo.event.wheelDelta},"mousewheel"):(Ta=function(){return-Zo.event.detail},"MozMousePixelScroll");Zo.color=et,et.prototype.toString=function(){return this.rgb()+""},Zo.hsl=rt;var Da=rt.prototype=new et;Da.brighter=function(n){return n=Math.pow(.7,arguments.length?n:1),new rt(this.h,this.s,this.l/n)},Da.darker=function(n){return n=Math.pow(.7,arguments.length?n:1),new rt(this.h,this.s,n*this.l)},Da.rgb=function(){return ut(this.h,this.s,this.l)},Zo.hcl=it;var Pa=it.prototype=new et;Pa.brighter=function(n){return new it(this.h,this.c,Math.min(100,this.l+Ua*(arguments.length?n:1)))},Pa.darker=function(n){return new it(this.h,this.c,Math.max(0,this.l-Ua*(arguments.length?n:1)))},Pa.rgb=function(){return ot(this.h,this.c,this.l).rgb()},Zo.lab=at;var Ua=18,ja=.95047,Ha=1,Fa=1.08883,Oa=at.prototype=new et;Oa.brighter=function(n){return new at(Math.min(100,this.l+Ua*(arguments.length?n:1)),this.a,this.b)},Oa.darker=function(n){return new at(Math.max(0,this.l-Ua*(arguments.length?n:1)),this.a,this.b)},Oa.rgb=function(){return ct(this.l,this.a,this.b)},Zo.rgb=gt;var Ya=gt.prototype=new et;Ya.brighter=function(n){n=Math.pow(.7,arguments.length?n:1);var t=this.r,e=this.g,r=this.b,u=30;return t||e||r?(t&&u>t&&(t=u),e&&u>e&&(e=u),r&&u>r&&(r=u),new gt(Math.min(255,t/n),Math.min(255,e/n),Math.min(255,r/n))):new gt(u,u,u)},Ya.darker=function(n){return n=Math.pow(.7,arguments.length?n:1),new gt(n*this.r,n*this.g,n*this.b)},Ya.hsl=function(){return yt(this.r,this.g,this.b)},Ya.toString=function(){return"#"+dt(this.r)+dt(this.g)+dt(this.b)};var Ia=Zo.map({aliceblue:15792383,antiquewhite:16444375,aqua:65535,aquamarine:8388564,azure:15794175,beige:16119260,bisque:16770244,black:0,blanchedalmond:16772045,blue:255,blueviolet:9055202,brown:10824234,burlywood:14596231,cadetblue:6266528,chartreuse:8388352,chocolate:13789470,coral:16744272,cornflowerblue:6591981,cornsilk:16775388,crimson:14423100,cyan:65535,darkblue:139,darkcyan:35723,darkgoldenrod:12092939,darkgray:11119017,darkgreen:25600,darkgrey:11119017,darkkhaki:12433259,darkmagenta:9109643,darkolivegreen:5597999,darkorange:16747520,darkorchid:10040012,darkred:9109504,darksalmon:15308410,darkseagreen:9419919,darkslateblue:4734347,darkslategray:3100495,darkslategrey:3100495,darkturquoise:52945,darkviolet:9699539,deeppink:16716947,deepskyblue:49151,dimgray:6908265,dimgrey:6908265,dodgerblue:2003199,firebrick:11674146,floralwhite:16775920,forestgreen:2263842,fuchsia:16711935,gainsboro:14474460,ghostwhite:16316671,gold:16766720,goldenrod:14329120,gray:8421504,green:32768,greenyellow:11403055,grey:8421504,honeydew:15794160,hotpink:16738740,indianred:13458524,indigo:4915330,ivory:16777200,khaki:15787660,lavender:15132410,lavenderblush:16773365,lawngreen:8190976,lemonchiffon:16775885,lightblue:11393254,lightcoral:15761536,lightcyan:14745599,lightgoldenrodyellow:16448210,lightgray:13882323,lightgreen:9498256,lightgrey:13882323,lightpink:16758465,lightsalmon:16752762,lightseagreen:2142890,lightskyblue:8900346,lightslategray:7833753,lightslategrey:7833753,lightsteelblue:11584734,lightyellow:16777184,lime:65280,limegreen:3329330,linen:16445670,magenta:16711935,maroon:8388608,mediumaquamarine:6737322,mediumblue:205,mediumorchid:12211667,mediumpurple:9662683,mediumseagreen:3978097,mediumslateblue:8087790,mediumspringgreen:64154,mediumturquoise:4772300,mediumvioletred:13047173,midnightblue:1644912,mintcream:16121850,mistyrose:16770273,moccasin:16770229,navajowhite:16768685,navy:128,oldlace:16643558,olive:8421376,olivedrab:7048739,orange:16753920,orangered:16729344,orchid:14315734,palegoldenrod:15657130,palegreen:10025880,paleturquoise:11529966,palevioletred:14381203,papayawhip:16773077,peachpuff:16767673,peru:13468991,pink:16761035,plum:14524637,powderblue:11591910,purple:8388736,red:16711680,rosybrown:12357519,royalblue:4286945,saddlebrown:9127187,salmon:16416882,sandybrown:16032864,seagreen:3050327,seashell:16774638,sienna:10506797,silver:12632256,skyblue:8900331,slateblue:6970061,slategray:7372944,slategrey:7372944,snow:16775930,springgreen:65407,steelblue:4620980,tan:13808780,teal:32896,thistle:14204888,tomato:16737095,turquoise:4251856,violet:15631086,wheat:16113331,white:16777215,whitesmoke:16119285,yellow:16776960,yellowgreen:10145074});Ia.forEach(function(n,t){Ia.set(n,pt(t))}),Zo.functor=bt,Zo.xhr=St(wt),Zo.dsv=function(n,t){function e(n,e,i){arguments.length<3&&(i=e,e=null);var o=kt(n,t,null==e?r:u(e),i);return o.row=function(n){return arguments.length?o.response(null==(e=n)?r:u(n)):e},o}function r(n){return e.parse(n.responseText)}function u(n){return function(t){return e.parse(t.responseText,n)}}function i(t){return t.map(o).join(n)}function o(n){return a.test(n)?'"'+n.replace(/\"/g,'""')+'"':n}var a=new RegExp('["'+n+"\n]"),c=n.charCodeAt(0);return e.parse=function(n,t){var r;return e.parseRows(n,function(n,e){if(r)return r(n,e-1);var u=new Function("d","return {"+n.map(function(n,t){return JSON.stringify(n)+": d["+t+"]"}).join(",")+"}");r=t?function(n,e){return t(u(n),e)}:u})},e.parseRows=function(n,t){function e(){if(l>=s)return o;if(u)return u=!1,i;var t=l;if(34===n.charCodeAt(t)){for(var e=t;e++<s;)if(34===n.charCodeAt(e)){if(34!==n.charCodeAt(e+1))break;++e}l=e+2;var r=n.charCodeAt(e+1);return 13===r?(u=!0,10===n.charCodeAt(e+2)&&++l):10===r&&(u=!0),n.substring(t+1,e).replace(/""/g,'"')}for(;s>l;){var r=n.charCodeAt(l++),a=1;if(10===r)u=!0;else if(13===r)u=!0,10===n.charCodeAt(l)&&(++l,++a);else if(r!==c)continue;return n.substring(t,l-a)}return n.substring(t)}for(var r,u,i={},o={},a=[],s=n.length,l=0,f=0;(r=e())!==o;){for(var h=[];r!==i&&r!==o;)h.push(r),r=e();(!t||(h=t(h,f++)))&&a.push(h)}return a},e.format=function(t){if(Array.isArray(t[0]))return e.formatRows(t);var r=new h,u=[];return t.forEach(function(n){for(var t in n)r.has(t)||u.push(r.add(t))}),[u.map(o).join(n)].concat(t.map(function(t){return u.map(function(n){return o(t[n])}).join(n)})).join("\n")},e.formatRows=function(n){return n.map(i).join("\n")},e},Zo.csv=Zo.dsv(",","text/csv"),Zo.tsv=Zo.dsv("	","text/tab-separated-values"),Zo.touch=function(n,t,e){if(arguments.length<3&&(e=t,t=x().changedTouches),t)for(var r,u=0,i=t.length;i>u;++u)if((r=t[u]).identifier===e)return Z(n,r)};var Za,Va,Xa,$a,Ba,Wa=Wo[p(Wo,"requestAnimationFrame")]||function(n){setTimeout(n,17)};Zo.timer=function(n,t,e){var r=arguments.length;2>r&&(t=0),3>r&&(e=Date.now());var u=e+t,i={c:n,t:u,f:!1,n:null};Va?Va.n=i:Za=i,Va=i,Xa||($a=clearTimeout($a),Xa=1,Wa(At))},Zo.timer.flush=function(){Ct(),Nt()},Zo.round=function(n,t){return t?Math.round(n*(t=Math.pow(10,t)))/t:Math.round(n)};var Ja=["y","z","a","f","p","n","\xb5","m","","k","M","G","T","P","E","Z","Y"].map(Lt);Zo.formatPrefix=function(n,t){var e=0;return n&&(0>n&&(n*=-1),t&&(n=Zo.round(n,zt(n,t))),e=1+Math.floor(1e-12+Math.log(n)/Math.LN10),e=Math.max(-24,Math.min(24,3*Math.floor((e-1)/3)))),Ja[8+e/3]};var Ga=/(?:([^{])?([<>=^]))?([+\- ])?([$#])?(0)?(\d+)?(,)?(\.-?\d+)?([a-z%])?/i,Ka=Zo.map({b:function(n){return n.toString(2)},c:function(n){return String.fromCharCode(n)},o:function(n){return n.toString(8)},x:function(n){return n.toString(16)},X:function(n){return n.toString(16).toUpperCase()},g:function(n,t){return n.toPrecision(t)},e:function(n,t){return n.toExponential(t)},f:function(n,t){return n.toFixed(t)},r:function(n,t){return(n=Zo.round(n,zt(n,t))).toFixed(Math.max(0,Math.min(20,zt(n*(1+1e-15),t))))}}),Qa=Zo.time={},nc=Date;Rt.prototype={getDate:function(){return this._.getUTCDate()},getDay:function(){return this._.getUTCDay()},getFullYear:function(){return this._.getUTCFullYear()},getHours:function(){return this._.getUTCHours()},getMilliseconds:function(){return this._.getUTCMilliseconds()},getMinutes:function(){return this._.getUTCMinutes()},getMonth:function(){return this._.getUTCMonth()},getSeconds:function(){return this._.getUTCSeconds()},getTime:function(){return this._.getTime()},getTimezoneOffset:function(){return 0},valueOf:function(){return this._.valueOf()},setDate:function(){tc.setUTCDate.apply(this._,arguments)},setDay:function(){tc.setUTCDay.apply(this._,arguments)},setFullYear:function(){tc.setUTCFullYear.apply(this._,arguments)},setHours:function(){tc.setUTCHours.apply(this._,arguments)},setMilliseconds:function(){tc.setUTCMilliseconds.apply(this._,arguments)},setMinutes:function(){tc.setUTCMinutes.apply(this._,arguments)},setMonth:function(){tc.setUTCMonth.apply(this._,arguments)},setSeconds:function(){tc.setUTCSeconds.apply(this._,arguments)},setTime:function(){tc.setTime.apply(this._,arguments)}};var tc=Date.prototype;Qa.year=Dt(function(n){return n=Qa.day(n),n.setMonth(0,1),n},function(n,t){n.setFullYear(n.getFullYear()+t)},function(n){return n.getFullYear()}),Qa.years=Qa.year.range,Qa.years.utc=Qa.year.utc.range,Qa.day=Dt(function(n){var t=new nc(2e3,0);return t.setFullYear(n.getFullYear(),n.getMonth(),n.getDate()),t},function(n,t){n.setDate(n.getDate()+t)},function(n){return n.getDate()-1}),Qa.days=Qa.day.range,Qa.days.utc=Qa.day.utc.range,Qa.dayOfYear=function(n){var t=Qa.year(n);return Math.floor((n-t-6e4*(n.getTimezoneOffset()-t.getTimezoneOffset()))/864e5)},["sunday","monday","tuesday","wednesday","thursday","friday","saturday"].forEach(function(n,t){t=7-t;var e=Qa[n]=Dt(function(n){return(n=Qa.day(n)).setDate(n.getDate()-(n.getDay()+t)%7),n},function(n,t){n.setDate(n.getDate()+7*Math.floor(t))},function(n){var e=Qa.year(n).getDay();return Math.floor((Qa.dayOfYear(n)+(e+t)%7)/7)-(e!==t)});Qa[n+"s"]=e.range,Qa[n+"s"].utc=e.utc.range,Qa[n+"OfYear"]=function(n){var e=Qa.year(n).getDay();return Math.floor((Qa.dayOfYear(n)+(e+t)%7)/7)}}),Qa.week=Qa.sunday,Qa.weeks=Qa.sunday.range,Qa.weeks.utc=Qa.sunday.utc.range,Qa.weekOfYear=Qa.sundayOfYear;var ec={"-":"",_:" ",0:"0"},rc=/^\s*\d+/,uc=/^%/;Zo.locale=function(n){return{numberFormat:Tt(n),timeFormat:Ut(n)}};var ic=Zo.locale({decimal:".",thousands:",",grouping:[3],currency:["$",""],dateTime:"%a %b %e %X %Y",date:"%m/%d/%Y",time:"%H:%M:%S",periods:["AM","PM"],days:["Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"],shortDays:["Sun","Mon","Tue","Wed","Thu","Fri","Sat"],months:["January","February","March","April","May","June","July","August","September","October","November","December"],shortMonths:["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]});Zo.format=ic.numberFormat,Zo.geo={},ue.prototype={s:0,t:0,add:function(n){ie(n,this.t,oc),ie(oc.s,this.s,this),this.s?this.t+=oc.t:this.s=oc.t},reset:function(){this.s=this.t=0},valueOf:function(){return this.s}};var oc=new ue;Zo.geo.stream=function(n,t){n&&ac.hasOwnProperty(n.type)?ac[n.type](n,t):oe(n,t)};var ac={Feature:function(n,t){oe(n.geometry,t)},FeatureCollection:function(n,t){for(var e=n.features,r=-1,u=e.length;++r<u;)oe(e[r].geometry,t)}},cc={Sphere:function(n,t){t.sphere()},Point:function(n,t){n=n.coordinates,t.point(n[0],n[1],n[2])},MultiPoint:function(n,t){for(var e=n.coordinates,r=-1,u=e.length;++r<u;)n=e[r],t.point(n[0],n[1],n[2])},LineString:function(n,t){ae(n.coordinates,t,0)},MultiLineString:function(n,t){for(var e=n.coordinates,r=-1,u=e.length;++r<u;)ae(e[r],t,0)},Polygon:function(n,t){ce(n.coordinates,t)},MultiPolygon:function(n,t){for(var e=n.coordinates,r=-1,u=e.length;++r<u;)ce(e[r],t)},GeometryCollection:function(n,t){for(var e=n.geometries,r=-1,u=e.length;++r<u;)oe(e[r],t)}};Zo.geo.area=function(n){return sc=0,Zo.geo.stream(n,fc),sc};var sc,lc=new ue,fc={sphere:function(){sc+=4*ba},point:v,lineStart:v,lineEnd:v,polygonStart:function(){lc.reset(),fc.lineStart=se},polygonEnd:function(){var n=2*lc;sc+=0>n?4*ba+n:n,fc.lineStart=fc.lineEnd=fc.point=v}};Zo.geo.bounds=function(){function n(n,t){x.push(M=[l=n,h=n]),f>t&&(f=t),t>g&&(g=t)}function t(t,e){var r=le([t*Aa,e*Aa]);if(m){var u=he(m,r),i=[u[1],-u[0],0],o=he(i,u);ve(o),o=de(o);var c=t-p,s=c>0?1:-1,v=o[0]*Ca*s,d=ua(c)>180;if(d^(v>s*p&&s*t>v)){var y=o[1]*Ca;y>g&&(g=y)}else if(v=(v+360)%360-180,d^(v>s*p&&s*t>v)){var y=-o[1]*Ca;f>y&&(f=y)}else f>e&&(f=e),e>g&&(g=e);d?p>t?a(l,t)>a(l,h)&&(h=t):a(t,h)>a(l,h)&&(l=t):h>=l?(l>t&&(l=t),t>h&&(h=t)):t>p?a(l,t)>a(l,h)&&(h=t):a(t,h)>a(l,h)&&(l=t)}else n(t,e);m=r,p=t}function e(){_.point=t}function r(){M[0]=l,M[1]=h,_.point=n,m=null}function u(n,e){if(m){var r=n-p;y+=ua(r)>180?r+(r>0?360:-360):r}else v=n,d=e;fc.point(n,e),t(n,e)}function i(){fc.lineStart()}function o(){u(v,d),fc.lineEnd(),ua(y)>ka&&(l=-(h=180)),M[0]=l,M[1]=h,m=null}function a(n,t){return(t-=n)<0?t+360:t}function c(n,t){return n[0]-t[0]}function s(n,t){return t[0]<=t[1]?t[0]<=n&&n<=t[1]:n<t[0]||t[1]<n}var l,f,h,g,p,v,d,m,y,x,M,_={point:n,lineStart:e,lineEnd:r,polygonStart:function(){_.point=u,_.lineStart=i,_.lineEnd=o,y=0,fc.polygonStart()},polygonEnd:function(){fc.polygonEnd(),_.point=n,_.lineStart=e,_.lineEnd=r,0>lc?(l=-(h=180),f=-(g=90)):y>ka?g=90:-ka>y&&(f=-90),M[0]=l,M[1]=h}};return function(n){g=h=-(l=f=1/0),x=[],Zo.geo.stream(n,_);var t=x.length;if(t){x.sort(c);for(var e,r=1,u=x[0],i=[u];t>r;++r)e=x[r],s(e[0],u)||s(e[1],u)?(a(u[0],e[1])>a(u[0],u[1])&&(u[1]=e[1]),a(e[0],u[1])>a(u[0],u[1])&&(u[0]=e[0])):i.push(u=e);
for(var o,e,p=-1/0,t=i.length-1,r=0,u=i[t];t>=r;u=e,++r)e=i[r],(o=a(u[1],e[0]))>p&&(p=o,l=e[0],h=u[1])}return x=M=null,1/0===l||1/0===f?[[0/0,0/0],[0/0,0/0]]:[[l,f],[h,g]]}}(),Zo.geo.centroid=function(n){hc=gc=pc=vc=dc=mc=yc=xc=Mc=_c=bc=0,Zo.geo.stream(n,wc);var t=Mc,e=_c,r=bc,u=t*t+e*e+r*r;return Ea>u&&(t=mc,e=yc,r=xc,ka>gc&&(t=pc,e=vc,r=dc),u=t*t+e*e+r*r,Ea>u)?[0/0,0/0]:[Math.atan2(e,t)*Ca,G(r/Math.sqrt(u))*Ca]};var hc,gc,pc,vc,dc,mc,yc,xc,Mc,_c,bc,wc={sphere:v,point:ye,lineStart:Me,lineEnd:_e,polygonStart:function(){wc.lineStart=be},polygonEnd:function(){wc.lineStart=Me}},Sc=Ae(we,Te,Re,[-ba,-ba/2]),kc=1e9;Zo.geo.clipExtent=function(){var n,t,e,r,u,i,o={stream:function(n){return u&&(u.valid=!1),u=i(n),u.valid=!0,u},extent:function(a){return arguments.length?(i=Ue(n=+a[0][0],t=+a[0][1],e=+a[1][0],r=+a[1][1]),u&&(u.valid=!1,u=null),o):[[n,t],[e,r]]}};return o.extent([[0,0],[960,500]])},(Zo.geo.conicEqualArea=function(){return He(Fe)}).raw=Fe,Zo.geo.albers=function(){return Zo.geo.conicEqualArea().rotate([96,0]).center([-.6,38.7]).parallels([29.5,45.5]).scale(1070)},Zo.geo.albersUsa=function(){function n(n){var i=n[0],o=n[1];return t=null,e(i,o),t||(r(i,o),t)||u(i,o),t}var t,e,r,u,i=Zo.geo.albers(),o=Zo.geo.conicEqualArea().rotate([154,0]).center([-2,58.5]).parallels([55,65]),a=Zo.geo.conicEqualArea().rotate([157,0]).center([-3,19.9]).parallels([8,18]),c={point:function(n,e){t=[n,e]}};return n.invert=function(n){var t=i.scale(),e=i.translate(),r=(n[0]-e[0])/t,u=(n[1]-e[1])/t;return(u>=.12&&.234>u&&r>=-.425&&-.214>r?o:u>=.166&&.234>u&&r>=-.214&&-.115>r?a:i).invert(n)},n.stream=function(n){var t=i.stream(n),e=o.stream(n),r=a.stream(n);return{point:function(n,u){t.point(n,u),e.point(n,u),r.point(n,u)},sphere:function(){t.sphere(),e.sphere(),r.sphere()},lineStart:function(){t.lineStart(),e.lineStart(),r.lineStart()},lineEnd:function(){t.lineEnd(),e.lineEnd(),r.lineEnd()},polygonStart:function(){t.polygonStart(),e.polygonStart(),r.polygonStart()},polygonEnd:function(){t.polygonEnd(),e.polygonEnd(),r.polygonEnd()}}},n.precision=function(t){return arguments.length?(i.precision(t),o.precision(t),a.precision(t),n):i.precision()},n.scale=function(t){return arguments.length?(i.scale(t),o.scale(.35*t),a.scale(t),n.translate(i.translate())):i.scale()},n.translate=function(t){if(!arguments.length)return i.translate();var s=i.scale(),l=+t[0],f=+t[1];return e=i.translate(t).clipExtent([[l-.455*s,f-.238*s],[l+.455*s,f+.238*s]]).stream(c).point,r=o.translate([l-.307*s,f+.201*s]).clipExtent([[l-.425*s+ka,f+.12*s+ka],[l-.214*s-ka,f+.234*s-ka]]).stream(c).point,u=a.translate([l-.205*s,f+.212*s]).clipExtent([[l-.214*s+ka,f+.166*s+ka],[l-.115*s-ka,f+.234*s-ka]]).stream(c).point,n},n.scale(1070)};var Ec,Ac,Cc,Nc,zc,Lc,Tc={point:v,lineStart:v,lineEnd:v,polygonStart:function(){Ac=0,Tc.lineStart=Oe},polygonEnd:function(){Tc.lineStart=Tc.lineEnd=Tc.point=v,Ec+=ua(Ac/2)}},qc={point:Ye,lineStart:v,lineEnd:v,polygonStart:v,polygonEnd:v},Rc={point:Ve,lineStart:Xe,lineEnd:$e,polygonStart:function(){Rc.lineStart=Be},polygonEnd:function(){Rc.point=Ve,Rc.lineStart=Xe,Rc.lineEnd=$e}};Zo.geo.path=function(){function n(n){return n&&("function"==typeof a&&i.pointRadius(+a.apply(this,arguments)),o&&o.valid||(o=u(i)),Zo.geo.stream(n,o)),i.result()}function t(){return o=null,n}var e,r,u,i,o,a=4.5;return n.area=function(n){return Ec=0,Zo.geo.stream(n,u(Tc)),Ec},n.centroid=function(n){return pc=vc=dc=mc=yc=xc=Mc=_c=bc=0,Zo.geo.stream(n,u(Rc)),bc?[Mc/bc,_c/bc]:xc?[mc/xc,yc/xc]:dc?[pc/dc,vc/dc]:[0/0,0/0]},n.bounds=function(n){return zc=Lc=-(Cc=Nc=1/0),Zo.geo.stream(n,u(qc)),[[Cc,Nc],[zc,Lc]]},n.projection=function(n){return arguments.length?(u=(e=n)?n.stream||Ge(n):wt,t()):e},n.context=function(n){return arguments.length?(i=null==(r=n)?new Ie:new We(n),"function"!=typeof a&&i.pointRadius(a),t()):r},n.pointRadius=function(t){return arguments.length?(a="function"==typeof t?t:(i.pointRadius(+t),+t),n):a},n.projection(Zo.geo.albersUsa()).context(null)},Zo.geo.transform=function(n){return{stream:function(t){var e=new Ke(t);for(var r in n)e[r]=n[r];return e}}},Ke.prototype={point:function(n,t){this.stream.point(n,t)},sphere:function(){this.stream.sphere()},lineStart:function(){this.stream.lineStart()},lineEnd:function(){this.stream.lineEnd()},polygonStart:function(){this.stream.polygonStart()},polygonEnd:function(){this.stream.polygonEnd()}},Zo.geo.projection=nr,Zo.geo.projectionMutator=tr,(Zo.geo.equirectangular=function(){return nr(rr)}).raw=rr.invert=rr,Zo.geo.rotation=function(n){function t(t){return t=n(t[0]*Aa,t[1]*Aa),t[0]*=Ca,t[1]*=Ca,t}return n=ir(n[0]%360*Aa,n[1]*Aa,n.length>2?n[2]*Aa:0),t.invert=function(t){return t=n.invert(t[0]*Aa,t[1]*Aa),t[0]*=Ca,t[1]*=Ca,t},t},ur.invert=rr,Zo.geo.circle=function(){function n(){var n="function"==typeof r?r.apply(this,arguments):r,t=ir(-n[0]*Aa,-n[1]*Aa,0).invert,u=[];return e(null,null,1,{point:function(n,e){u.push(n=t(n,e)),n[0]*=Ca,n[1]*=Ca}}),{type:"Polygon",coordinates:[u]}}var t,e,r=[0,0],u=6;return n.origin=function(t){return arguments.length?(r=t,n):r},n.angle=function(r){return arguments.length?(e=sr((t=+r)*Aa,u*Aa),n):t},n.precision=function(r){return arguments.length?(e=sr(t*Aa,(u=+r)*Aa),n):u},n.angle(90)},Zo.geo.distance=function(n,t){var e,r=(t[0]-n[0])*Aa,u=n[1]*Aa,i=t[1]*Aa,o=Math.sin(r),a=Math.cos(r),c=Math.sin(u),s=Math.cos(u),l=Math.sin(i),f=Math.cos(i);return Math.atan2(Math.sqrt((e=f*o)*e+(e=s*l-c*f*a)*e),c*l+s*f*a)},Zo.geo.graticule=function(){function n(){return{type:"MultiLineString",coordinates:t()}}function t(){return Zo.range(Math.ceil(i/d)*d,u,d).map(h).concat(Zo.range(Math.ceil(s/m)*m,c,m).map(g)).concat(Zo.range(Math.ceil(r/p)*p,e,p).filter(function(n){return ua(n%d)>ka}).map(l)).concat(Zo.range(Math.ceil(a/v)*v,o,v).filter(function(n){return ua(n%m)>ka}).map(f))}var e,r,u,i,o,a,c,s,l,f,h,g,p=10,v=p,d=90,m=360,y=2.5;return n.lines=function(){return t().map(function(n){return{type:"LineString",coordinates:n}})},n.outline=function(){return{type:"Polygon",coordinates:[h(i).concat(g(c).slice(1),h(u).reverse().slice(1),g(s).reverse().slice(1))]}},n.extent=function(t){return arguments.length?n.majorExtent(t).minorExtent(t):n.minorExtent()},n.majorExtent=function(t){return arguments.length?(i=+t[0][0],u=+t[1][0],s=+t[0][1],c=+t[1][1],i>u&&(t=i,i=u,u=t),s>c&&(t=s,s=c,c=t),n.precision(y)):[[i,s],[u,c]]},n.minorExtent=function(t){return arguments.length?(r=+t[0][0],e=+t[1][0],a=+t[0][1],o=+t[1][1],r>e&&(t=r,r=e,e=t),a>o&&(t=a,a=o,o=t),n.precision(y)):[[r,a],[e,o]]},n.step=function(t){return arguments.length?n.majorStep(t).minorStep(t):n.minorStep()},n.majorStep=function(t){return arguments.length?(d=+t[0],m=+t[1],n):[d,m]},n.minorStep=function(t){return arguments.length?(p=+t[0],v=+t[1],n):[p,v]},n.precision=function(t){return arguments.length?(y=+t,l=fr(a,o,90),f=hr(r,e,y),h=fr(s,c,90),g=hr(i,u,y),n):y},n.majorExtent([[-180,-90+ka],[180,90-ka]]).minorExtent([[-180,-80-ka],[180,80+ka]])},Zo.geo.greatArc=function(){function n(){return{type:"LineString",coordinates:[t||r.apply(this,arguments),e||u.apply(this,arguments)]}}var t,e,r=gr,u=pr;return n.distance=function(){return Zo.geo.distance(t||r.apply(this,arguments),e||u.apply(this,arguments))},n.source=function(e){return arguments.length?(r=e,t="function"==typeof e?null:e,n):r},n.target=function(t){return arguments.length?(u=t,e="function"==typeof t?null:t,n):u},n.precision=function(){return arguments.length?n:0},n},Zo.geo.interpolate=function(n,t){return vr(n[0]*Aa,n[1]*Aa,t[0]*Aa,t[1]*Aa)},Zo.geo.length=function(n){return Dc=0,Zo.geo.stream(n,Pc),Dc};var Dc,Pc={sphere:v,point:v,lineStart:dr,lineEnd:v,polygonStart:v,polygonEnd:v},Uc=mr(function(n){return Math.sqrt(2/(1+n))},function(n){return 2*Math.asin(n/2)});(Zo.geo.azimuthalEqualArea=function(){return nr(Uc)}).raw=Uc;var jc=mr(function(n){var t=Math.acos(n);return t&&t/Math.sin(t)},wt);(Zo.geo.azimuthalEquidistant=function(){return nr(jc)}).raw=jc,(Zo.geo.conicConformal=function(){return He(yr)}).raw=yr,(Zo.geo.conicEquidistant=function(){return He(xr)}).raw=xr;var Hc=mr(function(n){return 1/n},Math.atan);(Zo.geo.gnomonic=function(){return nr(Hc)}).raw=Hc,Mr.invert=function(n,t){return[n,2*Math.atan(Math.exp(t))-Sa]},(Zo.geo.mercator=function(){return _r(Mr)}).raw=Mr;var Fc=mr(function(){return 1},Math.asin);(Zo.geo.orthographic=function(){return nr(Fc)}).raw=Fc;var Oc=mr(function(n){return 1/(1+n)},function(n){return 2*Math.atan(n)});(Zo.geo.stereographic=function(){return nr(Oc)}).raw=Oc,br.invert=function(n,t){return[-t,2*Math.atan(Math.exp(n))-Sa]},(Zo.geo.transverseMercator=function(){var n=_r(br),t=n.center,e=n.rotate;return n.center=function(n){return n?t([-n[1],n[0]]):(n=t(),[n[1],-n[0]])},n.rotate=function(n){return n?e([n[0],n[1],n.length>2?n[2]+90:90]):(n=e(),[n[0],n[1],n[2]-90])},e([0,0,90])}).raw=br,Zo.geom={},Zo.geom.hull=function(n){function t(n){if(n.length<3)return[];var t,u=bt(e),i=bt(r),o=n.length,a=[],c=[];for(t=0;o>t;t++)a.push([+u.call(this,n[t],t),+i.call(this,n[t],t),t]);for(a.sort(Er),t=0;o>t;t++)c.push([a[t][0],-a[t][1]]);var s=kr(a),l=kr(c),f=l[0]===s[0],h=l[l.length-1]===s[s.length-1],g=[];for(t=s.length-1;t>=0;--t)g.push(n[a[s[t]][2]]);for(t=+f;t<l.length-h;++t)g.push(n[a[l[t]][2]]);return g}var e=wr,r=Sr;return arguments.length?t(n):(t.x=function(n){return arguments.length?(e=n,t):e},t.y=function(n){return arguments.length?(r=n,t):r},t)},Zo.geom.polygon=function(n){return sa(n,Yc),n};var Yc=Zo.geom.polygon.prototype=[];Yc.area=function(){for(var n,t=-1,e=this.length,r=this[e-1],u=0;++t<e;)n=r,r=this[t],u+=n[1]*r[0]-n[0]*r[1];return.5*u},Yc.centroid=function(n){var t,e,r=-1,u=this.length,i=0,o=0,a=this[u-1];for(arguments.length||(n=-1/(6*this.area()));++r<u;)t=a,a=this[r],e=t[0]*a[1]-a[0]*t[1],i+=(t[0]+a[0])*e,o+=(t[1]+a[1])*e;return[i*n,o*n]},Yc.clip=function(n){for(var t,e,r,u,i,o,a=Nr(n),c=-1,s=this.length-Nr(this),l=this[s-1];++c<s;){for(t=n.slice(),n.length=0,u=this[c],i=t[(r=t.length-a)-1],e=-1;++e<r;)o=t[e],Ar(o,l,u)?(Ar(i,l,u)||n.push(Cr(i,o,l,u)),n.push(o)):Ar(i,l,u)&&n.push(Cr(i,o,l,u)),i=o;a&&n.push(n[0]),l=u}return n};var Ic,Zc,Vc,Xc,$c,Bc=[],Wc=[];Ur.prototype.prepare=function(){for(var n,t=this.edges,e=t.length;e--;)n=t[e].edge,n.b&&n.a||t.splice(e,1);return t.sort(Hr),t.length},Wr.prototype={start:function(){return this.edge.l===this.site?this.edge.a:this.edge.b},end:function(){return this.edge.l===this.site?this.edge.b:this.edge.a}},Jr.prototype={insert:function(n,t){var e,r,u;if(n){if(t.P=n,t.N=n.N,n.N&&(n.N.P=t),n.N=t,n.R){for(n=n.R;n.L;)n=n.L;n.L=t}else n.R=t;e=n}else this._?(n=nu(this._),t.P=null,t.N=n,n.P=n.L=t,e=n):(t.P=t.N=null,this._=t,e=null);for(t.L=t.R=null,t.U=e,t.C=!0,n=t;e&&e.C;)r=e.U,e===r.L?(u=r.R,u&&u.C?(e.C=u.C=!1,r.C=!0,n=r):(n===e.R&&(Kr(this,e),n=e,e=n.U),e.C=!1,r.C=!0,Qr(this,r))):(u=r.L,u&&u.C?(e.C=u.C=!1,r.C=!0,n=r):(n===e.L&&(Qr(this,e),n=e,e=n.U),e.C=!1,r.C=!0,Kr(this,r))),e=n.U;this._.C=!1},remove:function(n){n.N&&(n.N.P=n.P),n.P&&(n.P.N=n.N),n.N=n.P=null;var t,e,r,u=n.U,i=n.L,o=n.R;if(e=i?o?nu(o):i:o,u?u.L===n?u.L=e:u.R=e:this._=e,i&&o?(r=e.C,e.C=n.C,e.L=i,i.U=e,e!==o?(u=e.U,e.U=n.U,n=e.R,u.L=n,e.R=o,o.U=e):(e.U=u,u=e,n=e.R)):(r=n.C,n=e),n&&(n.U=u),!r){if(n&&n.C)return n.C=!1,void 0;do{if(n===this._)break;if(n===u.L){if(t=u.R,t.C&&(t.C=!1,u.C=!0,Kr(this,u),t=u.R),t.L&&t.L.C||t.R&&t.R.C){t.R&&t.R.C||(t.L.C=!1,t.C=!0,Qr(this,t),t=u.R),t.C=u.C,u.C=t.R.C=!1,Kr(this,u),n=this._;break}}else if(t=u.L,t.C&&(t.C=!1,u.C=!0,Qr(this,u),t=u.L),t.L&&t.L.C||t.R&&t.R.C){t.L&&t.L.C||(t.R.C=!1,t.C=!0,Kr(this,t),t=u.L),t.C=u.C,u.C=t.L.C=!1,Qr(this,u),n=this._;break}t.C=!0,n=u,u=u.U}while(!n.C);n&&(n.C=!1)}}},Zo.geom.voronoi=function(n){function t(n){var t=new Array(n.length),r=a[0][0],u=a[0][1],i=a[1][0],o=a[1][1];return tu(e(n),a).cells.forEach(function(e,a){var c=e.edges,s=e.site,l=t[a]=c.length?c.map(function(n){var t=n.start();return[t.x,t.y]}):s.x>=r&&s.x<=i&&s.y>=u&&s.y<=o?[[r,o],[i,o],[i,u],[r,u]]:[];l.point=n[a]}),t}function e(n){return n.map(function(n,t){return{x:Math.round(i(n,t)/ka)*ka,y:Math.round(o(n,t)/ka)*ka,i:t}})}var r=wr,u=Sr,i=r,o=u,a=Jc;return n?t(n):(t.links=function(n){return tu(e(n)).edges.filter(function(n){return n.l&&n.r}).map(function(t){return{source:n[t.l.i],target:n[t.r.i]}})},t.triangles=function(n){var t=[];return tu(e(n)).cells.forEach(function(e,r){for(var u,i,o=e.site,a=e.edges.sort(Hr),c=-1,s=a.length,l=a[s-1].edge,f=l.l===o?l.r:l.l;++c<s;)u=l,i=f,l=a[c].edge,f=l.l===o?l.r:l.l,r<i.i&&r<f.i&&ru(o,i,f)<0&&t.push([n[r],n[i.i],n[f.i]])}),t},t.x=function(n){return arguments.length?(i=bt(r=n),t):r},t.y=function(n){return arguments.length?(o=bt(u=n),t):u},t.clipExtent=function(n){return arguments.length?(a=null==n?Jc:n,t):a===Jc?null:a},t.size=function(n){return arguments.length?t.clipExtent(n&&[[0,0],n]):a===Jc?null:a&&a[1]},t)};var Jc=[[-1e6,-1e6],[1e6,1e6]];Zo.geom.delaunay=function(n){return Zo.geom.voronoi().triangles(n)},Zo.geom.quadtree=function(n,t,e,r,u){function i(n){function i(n,t,e,r,u,i,o,a){if(!isNaN(e)&&!isNaN(r))if(n.leaf){var c=n.x,l=n.y;if(null!=c)if(ua(c-e)+ua(l-r)<.01)s(n,t,e,r,u,i,o,a);else{var f=n.point;n.x=n.y=n.point=null,s(n,f,c,l,u,i,o,a),s(n,t,e,r,u,i,o,a)}else n.x=e,n.y=r,n.point=t}else s(n,t,e,r,u,i,o,a)}function s(n,t,e,r,u,o,a,c){var s=.5*(u+a),l=.5*(o+c),f=e>=s,h=r>=l,g=(h<<1)+f;n.leaf=!1,n=n.nodes[g]||(n.nodes[g]=ou()),f?u=s:a=s,h?o=l:c=l,i(n,t,e,r,u,o,a,c)}var l,f,h,g,p,v,d,m,y,x=bt(a),M=bt(c);if(null!=t)v=t,d=e,m=r,y=u;else if(m=y=-(v=d=1/0),f=[],h=[],p=n.length,o)for(g=0;p>g;++g)l=n[g],l.x<v&&(v=l.x),l.y<d&&(d=l.y),l.x>m&&(m=l.x),l.y>y&&(y=l.y),f.push(l.x),h.push(l.y);else for(g=0;p>g;++g){var _=+x(l=n[g],g),b=+M(l,g);v>_&&(v=_),d>b&&(d=b),_>m&&(m=_),b>y&&(y=b),f.push(_),h.push(b)}var w=m-v,S=y-d;w>S?y=d+w:m=v+S;var k=ou();if(k.add=function(n){i(k,n,+x(n,++g),+M(n,g),v,d,m,y)},k.visit=function(n){au(n,k,v,d,m,y)},g=-1,null==t){for(;++g<p;)i(k,n[g],f[g],h[g],v,d,m,y);--g}else n.forEach(k.add);return f=h=n=l=null,k}var o,a=wr,c=Sr;return(o=arguments.length)?(a=uu,c=iu,3===o&&(u=e,r=t,e=t=0),i(n)):(i.x=function(n){return arguments.length?(a=n,i):a},i.y=function(n){return arguments.length?(c=n,i):c},i.extent=function(n){return arguments.length?(null==n?t=e=r=u=null:(t=+n[0][0],e=+n[0][1],r=+n[1][0],u=+n[1][1]),i):null==t?null:[[t,e],[r,u]]},i.size=function(n){return arguments.length?(null==n?t=e=r=u=null:(t=e=0,r=+n[0],u=+n[1]),i):null==t?null:[r-t,u-e]},i)},Zo.interpolateRgb=cu,Zo.interpolateObject=su,Zo.interpolateNumber=lu,Zo.interpolateString=fu;var Gc=/[-+]?(?:\d+\.?\d*|\.?\d+)(?:[eE][-+]?\d+)?/g,Kc=new RegExp(Gc.source,"g");Zo.interpolate=hu,Zo.interpolators=[function(n,t){var e=typeof t;return("string"===e?Ia.has(t)||/^(#|rgb\(|hsl\()/.test(t)?cu:fu:t instanceof et?cu:Array.isArray(t)?gu:"object"===e&&isNaN(t)?su:lu)(n,t)}],Zo.interpolateArray=gu;var Qc=function(){return wt},ns=Zo.map({linear:Qc,poly:Mu,quad:function(){return mu},cubic:function(){return yu},sin:function(){return _u},exp:function(){return bu},circle:function(){return wu},elastic:Su,back:ku,bounce:function(){return Eu}}),ts=Zo.map({"in":wt,out:vu,"in-out":du,"out-in":function(n){return du(vu(n))}});Zo.ease=function(n){var t=n.indexOf("-"),e=t>=0?n.substring(0,t):n,r=t>=0?n.substring(t+1):"in";return e=ns.get(e)||Qc,r=ts.get(r)||wt,pu(r(e.apply(null,Vo.call(arguments,1))))},Zo.interpolateHcl=Au,Zo.interpolateHsl=Cu,Zo.interpolateLab=Nu,Zo.interpolateRound=zu,Zo.transform=function(n){var t=$o.createElementNS(Zo.ns.prefix.svg,"g");return(Zo.transform=function(n){if(null!=n){t.setAttribute("transform",n);var e=t.transform.baseVal.consolidate()}return new Lu(e?e.matrix:es)})(n)},Lu.prototype.toString=function(){return"translate("+this.translate+")rotate("+this.rotate+")skewX("+this.skew+")scale("+this.scale+")"};var es={a:1,b:0,c:0,d:1,e:0,f:0};Zo.interpolateTransform=Du,Zo.layout={},Zo.layout.bundle=function(){return function(n){for(var t=[],e=-1,r=n.length;++e<r;)t.push(ju(n[e]));return t}},Zo.layout.chord=function(){function n(){var n,s,f,h,g,p={},v=[],d=Zo.range(i),m=[];for(e=[],r=[],n=0,h=-1;++h<i;){for(s=0,g=-1;++g<i;)s+=u[h][g];v.push(s),m.push(Zo.range(i)),n+=s}for(o&&d.sort(function(n,t){return o(v[n],v[t])}),a&&m.forEach(function(n,t){n.sort(function(n,e){return a(u[t][n],u[t][e])})}),n=(wa-l*i)/n,s=0,h=-1;++h<i;){for(f=s,g=-1;++g<i;){var y=d[h],x=m[y][g],M=u[y][x],_=s,b=s+=M*n;p[y+"-"+x]={index:y,subindex:x,startAngle:_,endAngle:b,value:M}}r[y]={index:y,startAngle:f,endAngle:s,value:(s-f)/n},s+=l}for(h=-1;++h<i;)for(g=h-1;++g<i;){var w=p[h+"-"+g],S=p[g+"-"+h];(w.value||S.value)&&e.push(w.value<S.value?{source:S,target:w}:{source:w,target:S})}c&&t()}function t(){e.sort(function(n,t){return c((n.source.value+n.target.value)/2,(t.source.value+t.target.value)/2)})}var e,r,u,i,o,a,c,s={},l=0;return s.matrix=function(n){return arguments.length?(i=(u=n)&&u.length,e=r=null,s):u},s.padding=function(n){return arguments.length?(l=n,e=r=null,s):l},s.sortGroups=function(n){return arguments.length?(o=n,e=r=null,s):o},s.sortSubgroups=function(n){return arguments.length?(a=n,e=null,s):a},s.sortChords=function(n){return arguments.length?(c=n,e&&t(),s):c},s.chords=function(){return e||n(),e},s.groups=function(){return r||n(),r},s},Zo.layout.force=function(){function n(n){return function(t,e,r,u){if(t.point!==n){var i=t.cx-n.x,o=t.cy-n.y,a=u-e,c=i*i+o*o;if(c>a*a/d){if(p>c){var s=t.charge/c;n.px-=i*s,n.py-=o*s}return!0}if(t.point&&c&&p>c){var s=t.pointCharge/c;n.px-=i*s,n.py-=o*s}}return!t.charge}}function t(n){n.px=Zo.event.x,n.py=Zo.event.y,a.resume()}var e,r,u,i,o,a={},c=Zo.dispatch("start","tick","end"),s=[1,1],l=.9,f=rs,h=us,g=-30,p=is,v=.1,d=.64,m=[],y=[];return a.tick=function(){if((r*=.99)<.005)return c.end({type:"end",alpha:r=0}),!0;var t,e,a,f,h,p,d,x,M,_=m.length,b=y.length;for(e=0;b>e;++e)a=y[e],f=a.source,h=a.target,x=h.x-f.x,M=h.y-f.y,(p=x*x+M*M)&&(p=r*i[e]*((p=Math.sqrt(p))-u[e])/p,x*=p,M*=p,h.x-=x*(d=f.weight/(h.weight+f.weight)),h.y-=M*d,f.x+=x*(d=1-d),f.y+=M*d);if((d=r*v)&&(x=s[0]/2,M=s[1]/2,e=-1,d))for(;++e<_;)a=m[e],a.x+=(x-a.x)*d,a.y+=(M-a.y)*d;if(g)for(Vu(t=Zo.geom.quadtree(m),r,o),e=-1;++e<_;)(a=m[e]).fixed||t.visit(n(a));for(e=-1;++e<_;)a=m[e],a.fixed?(a.x=a.px,a.y=a.py):(a.x-=(a.px-(a.px=a.x))*l,a.y-=(a.py-(a.py=a.y))*l);c.tick({type:"tick",alpha:r})},a.nodes=function(n){return arguments.length?(m=n,a):m},a.links=function(n){return arguments.length?(y=n,a):y},a.size=function(n){return arguments.length?(s=n,a):s},a.linkDistance=function(n){return arguments.length?(f="function"==typeof n?n:+n,a):f},a.distance=a.linkDistance,a.linkStrength=function(n){return arguments.length?(h="function"==typeof n?n:+n,a):h},a.friction=function(n){return arguments.length?(l=+n,a):l},a.charge=function(n){return arguments.length?(g="function"==typeof n?n:+n,a):g},a.chargeDistance=function(n){return arguments.length?(p=n*n,a):Math.sqrt(p)},a.gravity=function(n){return arguments.length?(v=+n,a):v},a.theta=function(n){return arguments.length?(d=n*n,a):Math.sqrt(d)},a.alpha=function(n){return arguments.length?(n=+n,r?r=n>0?n:0:n>0&&(c.start({type:"start",alpha:r=n}),Zo.timer(a.tick)),a):r},a.start=function(){function n(n,r){if(!e){for(e=new Array(c),a=0;c>a;++a)e[a]=[];for(a=0;s>a;++a){var u=y[a];e[u.source.index].push(u.target),e[u.target.index].push(u.source)}}for(var i,o=e[t],a=-1,s=o.length;++a<s;)if(!isNaN(i=o[a][n]))return i;return Math.random()*r}var t,e,r,c=m.length,l=y.length,p=s[0],v=s[1];for(t=0;c>t;++t)(r=m[t]).index=t,r.weight=0;for(t=0;l>t;++t)r=y[t],"number"==typeof r.source&&(r.source=m[r.source]),"number"==typeof r.target&&(r.target=m[r.target]),++r.source.weight,++r.target.weight;for(t=0;c>t;++t)r=m[t],isNaN(r.x)&&(r.x=n("x",p)),isNaN(r.y)&&(r.y=n("y",v)),isNaN(r.px)&&(r.px=r.x),isNaN(r.py)&&(r.py=r.y);if(u=[],"function"==typeof f)for(t=0;l>t;++t)u[t]=+f.call(this,y[t],t);else for(t=0;l>t;++t)u[t]=f;if(i=[],"function"==typeof h)for(t=0;l>t;++t)i[t]=+h.call(this,y[t],t);else for(t=0;l>t;++t)i[t]=h;if(o=[],"function"==typeof g)for(t=0;c>t;++t)o[t]=+g.call(this,m[t],t);else for(t=0;c>t;++t)o[t]=g;return a.resume()},a.resume=function(){return a.alpha(.1)},a.stop=function(){return a.alpha(0)},a.drag=function(){return e||(e=Zo.behavior.drag().origin(wt).on("dragstart.force",Ou).on("drag.force",t).on("dragend.force",Yu)),arguments.length?(this.on("mouseover.force",Iu).on("mouseout.force",Zu).call(e),void 0):e},Zo.rebind(a,c,"on")};var rs=20,us=1,is=1/0;Zo.layout.hierarchy=function(){function n(u){var i,o=[u],a=[];for(u.depth=0;null!=(i=o.pop());)if(a.push(i),(s=e.call(n,i,i.depth))&&(c=s.length)){for(var c,s,l;--c>=0;)o.push(l=s[c]),l.parent=i,l.depth=i.depth+1;r&&(i.value=0),i.children=s}else r&&(i.value=+r.call(n,i,i.depth)||0),delete i.children;return Bu(u,function(n){var e,u;t&&(e=n.children)&&e.sort(t),r&&(u=n.parent)&&(u.value+=n.value)}),a}var t=Gu,e=Wu,r=Ju;return n.sort=function(e){return arguments.length?(t=e,n):t},n.children=function(t){return arguments.length?(e=t,n):e},n.value=function(t){return arguments.length?(r=t,n):r},n.revalue=function(t){return r&&($u(t,function(n){n.children&&(n.value=0)}),Bu(t,function(t){var e;t.children||(t.value=+r.call(n,t,t.depth)||0),(e=t.parent)&&(e.value+=t.value)})),t},n},Zo.layout.partition=function(){function n(t,e,r,u){var i=t.children;if(t.x=e,t.y=t.depth*u,t.dx=r,t.dy=u,i&&(o=i.length)){var o,a,c,s=-1;for(r=t.value?r/t.value:0;++s<o;)n(a=i[s],e,c=a.value*r,u),e+=c}}function t(n){var e=n.children,r=0;if(e&&(u=e.length))for(var u,i=-1;++i<u;)r=Math.max(r,t(e[i]));return 1+r}function e(e,i){var o=r.call(this,e,i);return n(o[0],0,u[0],u[1]/t(o[0])),o}var r=Zo.layout.hierarchy(),u=[1,1];return e.size=function(n){return arguments.length?(u=n,e):u},Xu(e,r)},Zo.layout.pie=function(){function n(i){var o=i.map(function(e,r){return+t.call(n,e,r)}),a=+("function"==typeof r?r.apply(this,arguments):r),c=(("function"==typeof u?u.apply(this,arguments):u)-a)/Zo.sum(o),s=Zo.range(i.length);null!=e&&s.sort(e===os?function(n,t){return o[t]-o[n]}:function(n,t){return e(i[n],i[t])});var l=[];return s.forEach(function(n){var t;l[n]={data:i[n],value:t=o[n],startAngle:a,endAngle:a+=t*c}}),l}var t=Number,e=os,r=0,u=wa;return n.value=function(e){return arguments.length?(t=e,n):t},n.sort=function(t){return arguments.length?(e=t,n):e},n.startAngle=function(t){return arguments.length?(r=t,n):r},n.endAngle=function(t){return arguments.length?(u=t,n):u},n};var os={};Zo.layout.stack=function(){function n(a,c){var s=a.map(function(e,r){return t.call(n,e,r)}),l=s.map(function(t){return t.map(function(t,e){return[i.call(n,t,e),o.call(n,t,e)]})}),f=e.call(n,l,c);s=Zo.permute(s,f),l=Zo.permute(l,f);var h,g,p,v=r.call(n,l,c),d=s.length,m=s[0].length;for(g=0;m>g;++g)for(u.call(n,s[0][g],p=v[g],l[0][g][1]),h=1;d>h;++h)u.call(n,s[h][g],p+=l[h-1][g][1],l[h][g][1]);return a}var t=wt,e=ei,r=ri,u=ti,i=Qu,o=ni;return n.values=function(e){return arguments.length?(t=e,n):t},n.order=function(t){return arguments.length?(e="function"==typeof t?t:as.get(t)||ei,n):e},n.offset=function(t){return arguments.length?(r="function"==typeof t?t:cs.get(t)||ri,n):r},n.x=function(t){return arguments.length?(i=t,n):i},n.y=function(t){return arguments.length?(o=t,n):o},n.out=function(t){return arguments.length?(u=t,n):u},n};var as=Zo.map({"inside-out":function(n){var t,e,r=n.length,u=n.map(ui),i=n.map(ii),o=Zo.range(r).sort(function(n,t){return u[n]-u[t]}),a=0,c=0,s=[],l=[];for(t=0;r>t;++t)e=o[t],c>a?(a+=i[e],s.push(e)):(c+=i[e],l.push(e));return l.reverse().concat(s)},reverse:function(n){return Zo.range(n.length).reverse()},"default":ei}),cs=Zo.map({silhouette:function(n){var t,e,r,u=n.length,i=n[0].length,o=[],a=0,c=[];for(e=0;i>e;++e){for(t=0,r=0;u>t;t++)r+=n[t][e][1];r>a&&(a=r),o.push(r)}for(e=0;i>e;++e)c[e]=(a-o[e])/2;return c},wiggle:function(n){var t,e,r,u,i,o,a,c,s,l=n.length,f=n[0],h=f.length,g=[];for(g[0]=c=s=0,e=1;h>e;++e){for(t=0,u=0;l>t;++t)u+=n[t][e][1];for(t=0,i=0,a=f[e][0]-f[e-1][0];l>t;++t){for(r=0,o=(n[t][e][1]-n[t][e-1][1])/(2*a);t>r;++r)o+=(n[r][e][1]-n[r][e-1][1])/a;i+=o*n[t][e][1]}g[e]=c-=u?i/u*a:0,s>c&&(s=c)}for(e=0;h>e;++e)g[e]-=s;return g},expand:function(n){var t,e,r,u=n.length,i=n[0].length,o=1/u,a=[];for(e=0;i>e;++e){for(t=0,r=0;u>t;t++)r+=n[t][e][1];if(r)for(t=0;u>t;t++)n[t][e][1]/=r;else for(t=0;u>t;t++)n[t][e][1]=o}for(e=0;i>e;++e)a[e]=0;return a},zero:ri});Zo.layout.histogram=function(){function n(n,i){for(var o,a,c=[],s=n.map(e,this),l=r.call(this,s,i),f=u.call(this,l,s,i),i=-1,h=s.length,g=f.length-1,p=t?1:1/h;++i<g;)o=c[i]=[],o.dx=f[i+1]-(o.x=f[i]),o.y=0;if(g>0)for(i=-1;++i<h;)a=s[i],a>=l[0]&&a<=l[1]&&(o=c[Zo.bisect(f,a,1,g)-1],o.y+=p,o.push(n[i]));return c}var t=!0,e=Number,r=si,u=ai;return n.value=function(t){return arguments.length?(e=t,n):e},n.range=function(t){return arguments.length?(r=bt(t),n):r},n.bins=function(t){return arguments.length?(u="number"==typeof t?function(n){return ci(n,t)}:bt(t),n):u},n.frequency=function(e){return arguments.length?(t=!!e,n):t},n},Zo.layout.pack=function(){function n(n,i){var o=e.call(this,n,i),a=o[0],c=u[0],s=u[1],l=null==t?Math.sqrt:"function"==typeof t?t:function(){return t};if(a.x=a.y=0,Bu(a,function(n){n.r=+l(n.value)}),Bu(a,pi),r){var f=r*(t?1:Math.max(2*a.r/c,2*a.r/s))/2;Bu(a,function(n){n.r+=f}),Bu(a,pi),Bu(a,function(n){n.r-=f})}return mi(a,c/2,s/2,t?1:1/Math.max(2*a.r/c,2*a.r/s)),o}var t,e=Zo.layout.hierarchy().sort(li),r=0,u=[1,1];return n.size=function(t){return arguments.length?(u=t,n):u},n.radius=function(e){return arguments.length?(t=null==e||"function"==typeof e?e:+e,n):t},n.padding=function(t){return arguments.length?(r=+t,n):r},Xu(n,e)},Zo.layout.tree=function(){function n(n,u){var l=o.call(this,n,u),f=l[0],h=t(f);if(Bu(h,e),h.parent.m=-h.z,$u(h,r),s)$u(f,i);else{var g=f,p=f,v=f;$u(f,function(n){n.x<g.x&&(g=n),n.x>p.x&&(p=n),n.depth>v.depth&&(v=n)});var d=a(g,p)/2-g.x,m=c[0]/(p.x+a(p,g)/2+d),y=c[1]/(v.depth||1);$u(f,function(n){n.x=(n.x+d)*m,n.y=n.depth*y})}return l}function t(n){for(var t,e={A:null,children:[n]},r=[e];null!=(t=r.pop());)for(var u,i=t.children,o=0,a=i.length;a>o;++o)r.push((i[o]=u={_:i[o],parent:t,children:(u=i[o].children)&&u.slice()||[],A:null,a:null,z:0,m:0,c:0,s:0,t:null,i:o}).a=u);return e.children[0]}function e(n){var t=n.children,e=n.parent.children,r=n.i?e[n.i-1]:null;if(t.length){wi(n);var i=(t[0].z+t[t.length-1].z)/2;r?(n.z=r.z+a(n._,r._),n.m=n.z-i):n.z=i}else r&&(n.z=r.z+a(n._,r._));n.parent.A=u(n,r,n.parent.A||e[0])}function r(n){n._.x=n.z+n.parent.m,n.m+=n.parent.m}function u(n,t,e){if(t){for(var r,u=n,i=n,o=t,c=u.parent.children[0],s=u.m,l=i.m,f=o.m,h=c.m;o=_i(o),u=Mi(u),o&&u;)c=Mi(c),i=_i(i),i.a=n,r=o.z+f-u.z-s+a(o._,u._),r>0&&(bi(Si(o,n,e),n,r),s+=r,l+=r),f+=o.m,s+=u.m,h+=c.m,l+=i.m;o&&!_i(i)&&(i.t=o,i.m+=f-l),u&&!Mi(c)&&(c.t=u,c.m+=s-h,e=n)}return e}function i(n){n.x*=c[0],n.y=n.depth*c[1]}var o=Zo.layout.hierarchy().sort(null).value(null),a=xi,c=[1,1],s=null;return n.separation=function(t){return arguments.length?(a=t,n):a},n.size=function(t){return arguments.length?(s=null==(c=t)?i:null,n):s?null:c},n.nodeSize=function(t){return arguments.length?(s=null==(c=t)?null:i,n):s?c:null},Xu(n,o)},Zo.layout.cluster=function(){function n(n,i){var o,a=t.call(this,n,i),c=a[0],s=0;Bu(c,function(n){var t=n.children;t&&t.length?(n.x=Ei(t),n.y=ki(t)):(n.x=o?s+=e(n,o):0,n.y=0,o=n)});var l=Ai(c),f=Ci(c),h=l.x-e(l,f)/2,g=f.x+e(f,l)/2;return Bu(c,u?function(n){n.x=(n.x-c.x)*r[0],n.y=(c.y-n.y)*r[1]}:function(n){n.x=(n.x-h)/(g-h)*r[0],n.y=(1-(c.y?n.y/c.y:1))*r[1]}),a}var t=Zo.layout.hierarchy().sort(null).value(null),e=xi,r=[1,1],u=!1;return n.separation=function(t){return arguments.length?(e=t,n):e},n.size=function(t){return arguments.length?(u=null==(r=t),n):u?null:r},n.nodeSize=function(t){return arguments.length?(u=null!=(r=t),n):u?r:null},Xu(n,t)},Zo.layout.treemap=function(){function n(n,t){for(var e,r,u=-1,i=n.length;++u<i;)r=(e=n[u]).value*(0>t?0:t),e.area=isNaN(r)||0>=r?0:r}function t(e){var i=e.children;if(i&&i.length){var o,a,c,s=f(e),l=[],h=i.slice(),p=1/0,v="slice"===g?s.dx:"dice"===g?s.dy:"slice-dice"===g?1&e.depth?s.dy:s.dx:Math.min(s.dx,s.dy);for(n(h,s.dx*s.dy/e.value),l.area=0;(c=h.length)>0;)l.push(o=h[c-1]),l.area+=o.area,"squarify"!==g||(a=r(l,v))<=p?(h.pop(),p=a):(l.area-=l.pop().area,u(l,v,s,!1),v=Math.min(s.dx,s.dy),l.length=l.area=0,p=1/0);l.length&&(u(l,v,s,!0),l.length=l.area=0),i.forEach(t)}}function e(t){var r=t.children;if(r&&r.length){var i,o=f(t),a=r.slice(),c=[];for(n(a,o.dx*o.dy/t.value),c.area=0;i=a.pop();)c.push(i),c.area+=i.area,null!=i.z&&(u(c,i.z?o.dx:o.dy,o,!a.length),c.length=c.area=0);r.forEach(e)}}function r(n,t){for(var e,r=n.area,u=0,i=1/0,o=-1,a=n.length;++o<a;)(e=n[o].area)&&(i>e&&(i=e),e>u&&(u=e));return r*=r,t*=t,r?Math.max(t*u*p/r,r/(t*i*p)):1/0}function u(n,t,e,r){var u,i=-1,o=n.length,a=e.x,s=e.y,l=t?c(n.area/t):0;if(t==e.dx){for((r||l>e.dy)&&(l=e.dy);++i<o;)u=n[i],u.x=a,u.y=s,u.dy=l,a+=u.dx=Math.min(e.x+e.dx-a,l?c(u.area/l):0);u.z=!0,u.dx+=e.x+e.dx-a,e.y+=l,e.dy-=l}else{for((r||l>e.dx)&&(l=e.dx);++i<o;)u=n[i],u.x=a,u.y=s,u.dx=l,s+=u.dy=Math.min(e.y+e.dy-s,l?c(u.area/l):0);u.z=!1,u.dy+=e.y+e.dy-s,e.x+=l,e.dx-=l}}function i(r){var u=o||a(r),i=u[0];return i.x=0,i.y=0,i.dx=s[0],i.dy=s[1],o&&a.revalue(i),n([i],i.dx*i.dy/i.value),(o?e:t)(i),h&&(o=u),u}var o,a=Zo.layout.hierarchy(),c=Math.round,s=[1,1],l=null,f=Ni,h=!1,g="squarify",p=.5*(1+Math.sqrt(5));return i.size=function(n){return arguments.length?(s=n,i):s},i.padding=function(n){function t(t){var e=n.call(i,t,t.depth);return null==e?Ni(t):zi(t,"number"==typeof e?[e,e,e,e]:e)}function e(t){return zi(t,n)}if(!arguments.length)return l;var r;return f=null==(l=n)?Ni:"function"==(r=typeof n)?t:"number"===r?(n=[n,n,n,n],e):e,i},i.round=function(n){return arguments.length?(c=n?Math.round:Number,i):c!=Number},i.sticky=function(n){return arguments.length?(h=n,o=null,i):h},i.ratio=function(n){return arguments.length?(p=n,i):p},i.mode=function(n){return arguments.length?(g=n+"",i):g},Xu(i,a)},Zo.random={normal:function(n,t){var e=arguments.length;return 2>e&&(t=1),1>e&&(n=0),function(){var e,r,u;do e=2*Math.random()-1,r=2*Math.random()-1,u=e*e+r*r;while(!u||u>1);return n+t*e*Math.sqrt(-2*Math.log(u)/u)}},logNormal:function(){var n=Zo.random.normal.apply(Zo,arguments);return function(){return Math.exp(n())}},bates:function(n){var t=Zo.random.irwinHall(n);return function(){return t()/n}},irwinHall:function(n){return function(){for(var t=0,e=0;n>e;e++)t+=Math.random();return t}}},Zo.scale={};var ss={floor:wt,ceil:wt};Zo.scale.linear=function(){return Ui([0,1],[0,1],hu,!1)};var ls={s:1,g:1,p:1,r:1,e:1};Zo.scale.log=function(){return Vi(Zo.scale.linear().domain([0,1]),10,!0,[1,10])};var fs=Zo.format(".0e"),hs={floor:function(n){return-Math.ceil(-n)},ceil:function(n){return-Math.floor(-n)}};Zo.scale.pow=function(){return Xi(Zo.scale.linear(),1,[0,1])},Zo.scale.sqrt=function(){return Zo.scale.pow().exponent(.5)},Zo.scale.ordinal=function(){return Bi([],{t:"range",a:[[]]})},Zo.scale.category10=function(){return Zo.scale.ordinal().range(gs)},Zo.scale.category20=function(){return Zo.scale.ordinal().range(ps)},Zo.scale.category20b=function(){return Zo.scale.ordinal().range(vs)},Zo.scale.category20c=function(){return Zo.scale.ordinal().range(ds)};var gs=[2062260,16744206,2924588,14034728,9725885,9197131,14907330,8355711,12369186,1556175].map(vt),ps=[2062260,11454440,16744206,16759672,2924588,10018698,14034728,16750742,9725885,12955861,9197131,12885140,14907330,16234194,8355711,13092807,12369186,14408589,1556175,10410725].map(vt),vs=[3750777,5395619,7040719,10264286,6519097,9216594,11915115,13556636,9202993,12426809,15186514,15190932,8666169,11356490,14049643,15177372,8077683,10834324,13528509,14589654].map(vt),ds=[3244733,7057110,10406625,13032431,15095053,16616764,16625259,16634018,3253076,7652470,10607003,13101504,7695281,10394312,12369372,14342891,6513507,9868950,12434877,14277081].map(vt);Zo.scale.quantile=function(){return Wi([],[])},Zo.scale.quantize=function(){return Ji(0,1,[0,1])},Zo.scale.threshold=function(){return Gi([.5],[0,1])},Zo.scale.identity=function(){return Ki([0,1])},Zo.svg={},Zo.svg.arc=function(){function n(){var n=t.apply(this,arguments),i=e.apply(this,arguments),o=r.apply(this,arguments)+ms,a=u.apply(this,arguments)+ms,c=(o>a&&(c=o,o=a,a=c),a-o),s=ba>c?"0":"1",l=Math.cos(o),f=Math.sin(o),h=Math.cos(a),g=Math.sin(a);
return c>=ys?n?"M0,"+i+"A"+i+","+i+" 0 1,1 0,"+-i+"A"+i+","+i+" 0 1,1 0,"+i+"M0,"+n+"A"+n+","+n+" 0 1,0 0,"+-n+"A"+n+","+n+" 0 1,0 0,"+n+"Z":"M0,"+i+"A"+i+","+i+" 0 1,1 0,"+-i+"A"+i+","+i+" 0 1,1 0,"+i+"Z":n?"M"+i*l+","+i*f+"A"+i+","+i+" 0 "+s+",1 "+i*h+","+i*g+"L"+n*h+","+n*g+"A"+n+","+n+" 0 "+s+",0 "+n*l+","+n*f+"Z":"M"+i*l+","+i*f+"A"+i+","+i+" 0 "+s+",1 "+i*h+","+i*g+"L0,0"+"Z"}var t=Qi,e=no,r=to,u=eo;return n.innerRadius=function(e){return arguments.length?(t=bt(e),n):t},n.outerRadius=function(t){return arguments.length?(e=bt(t),n):e},n.startAngle=function(t){return arguments.length?(r=bt(t),n):r},n.endAngle=function(t){return arguments.length?(u=bt(t),n):u},n.centroid=function(){var n=(t.apply(this,arguments)+e.apply(this,arguments))/2,i=(r.apply(this,arguments)+u.apply(this,arguments))/2+ms;return[Math.cos(i)*n,Math.sin(i)*n]},n};var ms=-Sa,ys=wa-ka;Zo.svg.line=function(){return ro(wt)};var xs=Zo.map({linear:uo,"linear-closed":io,step:oo,"step-before":ao,"step-after":co,basis:po,"basis-open":vo,"basis-closed":mo,bundle:yo,cardinal:fo,"cardinal-open":so,"cardinal-closed":lo,monotone:So});xs.forEach(function(n,t){t.key=n,t.closed=/-closed$/.test(n)});var Ms=[0,2/3,1/3,0],_s=[0,1/3,2/3,0],bs=[0,1/6,2/3,1/6];Zo.svg.line.radial=function(){var n=ro(ko);return n.radius=n.x,delete n.x,n.angle=n.y,delete n.y,n},ao.reverse=co,co.reverse=ao,Zo.svg.area=function(){return Eo(wt)},Zo.svg.area.radial=function(){var n=Eo(ko);return n.radius=n.x,delete n.x,n.innerRadius=n.x0,delete n.x0,n.outerRadius=n.x1,delete n.x1,n.angle=n.y,delete n.y,n.startAngle=n.y0,delete n.y0,n.endAngle=n.y1,delete n.y1,n},Zo.svg.chord=function(){function n(n,a){var c=t(this,i,n,a),s=t(this,o,n,a);return"M"+c.p0+r(c.r,c.p1,c.a1-c.a0)+(e(c,s)?u(c.r,c.p1,c.r,c.p0):u(c.r,c.p1,s.r,s.p0)+r(s.r,s.p1,s.a1-s.a0)+u(s.r,s.p1,c.r,c.p0))+"Z"}function t(n,t,e,r){var u=t.call(n,e,r),i=a.call(n,u,r),o=c.call(n,u,r)+ms,l=s.call(n,u,r)+ms;return{r:i,a0:o,a1:l,p0:[i*Math.cos(o),i*Math.sin(o)],p1:[i*Math.cos(l),i*Math.sin(l)]}}function e(n,t){return n.a0==t.a0&&n.a1==t.a1}function r(n,t,e){return"A"+n+","+n+" 0 "+ +(e>ba)+",1 "+t}function u(n,t,e,r){return"Q 0,0 "+r}var i=gr,o=pr,a=Ao,c=to,s=eo;return n.radius=function(t){return arguments.length?(a=bt(t),n):a},n.source=function(t){return arguments.length?(i=bt(t),n):i},n.target=function(t){return arguments.length?(o=bt(t),n):o},n.startAngle=function(t){return arguments.length?(c=bt(t),n):c},n.endAngle=function(t){return arguments.length?(s=bt(t),n):s},n},Zo.svg.diagonal=function(){function n(n,u){var i=t.call(this,n,u),o=e.call(this,n,u),a=(i.y+o.y)/2,c=[i,{x:i.x,y:a},{x:o.x,y:a},o];return c=c.map(r),"M"+c[0]+"C"+c[1]+" "+c[2]+" "+c[3]}var t=gr,e=pr,r=Co;return n.source=function(e){return arguments.length?(t=bt(e),n):t},n.target=function(t){return arguments.length?(e=bt(t),n):e},n.projection=function(t){return arguments.length?(r=t,n):r},n},Zo.svg.diagonal.radial=function(){var n=Zo.svg.diagonal(),t=Co,e=n.projection;return n.projection=function(n){return arguments.length?e(No(t=n)):t},n},Zo.svg.symbol=function(){function n(n,r){return(ws.get(t.call(this,n,r))||To)(e.call(this,n,r))}var t=Lo,e=zo;return n.type=function(e){return arguments.length?(t=bt(e),n):t},n.size=function(t){return arguments.length?(e=bt(t),n):e},n};var ws=Zo.map({circle:To,cross:function(n){var t=Math.sqrt(n/5)/2;return"M"+-3*t+","+-t+"H"+-t+"V"+-3*t+"H"+t+"V"+-t+"H"+3*t+"V"+t+"H"+t+"V"+3*t+"H"+-t+"V"+t+"H"+-3*t+"Z"},diamond:function(n){var t=Math.sqrt(n/(2*As)),e=t*As;return"M0,"+-t+"L"+e+",0"+" 0,"+t+" "+-e+",0"+"Z"},square:function(n){var t=Math.sqrt(n)/2;return"M"+-t+","+-t+"L"+t+","+-t+" "+t+","+t+" "+-t+","+t+"Z"},"triangle-down":function(n){var t=Math.sqrt(n/Es),e=t*Es/2;return"M0,"+e+"L"+t+","+-e+" "+-t+","+-e+"Z"},"triangle-up":function(n){var t=Math.sqrt(n/Es),e=t*Es/2;return"M0,"+-e+"L"+t+","+e+" "+-t+","+e+"Z"}});Zo.svg.symbolTypes=ws.keys();var Ss,ks,Es=Math.sqrt(3),As=Math.tan(30*Aa),Cs=[],Ns=0;Cs.call=pa.call,Cs.empty=pa.empty,Cs.node=pa.node,Cs.size=pa.size,Zo.transition=function(n){return arguments.length?Ss?n.transition():n:ma.transition()},Zo.transition.prototype=Cs,Cs.select=function(n){var t,e,r,u=this.id,i=[];n=b(n);for(var o=-1,a=this.length;++o<a;){i.push(t=[]);for(var c=this[o],s=-1,l=c.length;++s<l;)(r=c[s])&&(e=n.call(r,r.__data__,s,o))?("__data__"in r&&(e.__data__=r.__data__),Po(e,s,u,r.__transition__[u]),t.push(e)):t.push(null)}return qo(i,u)},Cs.selectAll=function(n){var t,e,r,u,i,o=this.id,a=[];n=w(n);for(var c=-1,s=this.length;++c<s;)for(var l=this[c],f=-1,h=l.length;++f<h;)if(r=l[f]){i=r.__transition__[o],e=n.call(r,r.__data__,f,c),a.push(t=[]);for(var g=-1,p=e.length;++g<p;)(u=e[g])&&Po(u,g,o,i),t.push(u)}return qo(a,o)},Cs.filter=function(n){var t,e,r,u=[];"function"!=typeof n&&(n=R(n));for(var i=0,o=this.length;o>i;i++){u.push(t=[]);for(var e=this[i],a=0,c=e.length;c>a;a++)(r=e[a])&&n.call(r,r.__data__,a,i)&&t.push(r)}return qo(u,this.id)},Cs.tween=function(n,t){var e=this.id;return arguments.length<2?this.node().__transition__[e].tween.get(n):P(this,null==t?function(t){t.__transition__[e].tween.remove(n)}:function(r){r.__transition__[e].tween.set(n,t)})},Cs.attr=function(n,t){function e(){this.removeAttribute(a)}function r(){this.removeAttributeNS(a.space,a.local)}function u(n){return null==n?e:(n+="",function(){var t,e=this.getAttribute(a);return e!==n&&(t=o(e,n),function(n){this.setAttribute(a,t(n))})})}function i(n){return null==n?r:(n+="",function(){var t,e=this.getAttributeNS(a.space,a.local);return e!==n&&(t=o(e,n),function(n){this.setAttributeNS(a.space,a.local,t(n))})})}if(arguments.length<2){for(t in n)this.attr(t,n[t]);return this}var o="transform"==n?Du:hu,a=Zo.ns.qualify(n);return Ro(this,"attr."+n,t,a.local?i:u)},Cs.attrTween=function(n,t){function e(n,e){var r=t.call(this,n,e,this.getAttribute(u));return r&&function(n){this.setAttribute(u,r(n))}}function r(n,e){var r=t.call(this,n,e,this.getAttributeNS(u.space,u.local));return r&&function(n){this.setAttributeNS(u.space,u.local,r(n))}}var u=Zo.ns.qualify(n);return this.tween("attr."+n,u.local?r:e)},Cs.style=function(n,t,e){function r(){this.style.removeProperty(n)}function u(t){return null==t?r:(t+="",function(){var r,u=Wo.getComputedStyle(this,null).getPropertyValue(n);return u!==t&&(r=hu(u,t),function(t){this.style.setProperty(n,r(t),e)})})}var i=arguments.length;if(3>i){if("string"!=typeof n){2>i&&(t="");for(e in n)this.style(e,n[e],t);return this}e=""}return Ro(this,"style."+n,t,u)},Cs.styleTween=function(n,t,e){function r(r,u){var i=t.call(this,r,u,Wo.getComputedStyle(this,null).getPropertyValue(n));return i&&function(t){this.style.setProperty(n,i(t),e)}}return arguments.length<3&&(e=""),this.tween("style."+n,r)},Cs.text=function(n){return Ro(this,"text",n,Do)},Cs.remove=function(){return this.each("end.transition",function(){var n;this.__transition__.count<2&&(n=this.parentNode)&&n.removeChild(this)})},Cs.ease=function(n){var t=this.id;return arguments.length<1?this.node().__transition__[t].ease:("function"!=typeof n&&(n=Zo.ease.apply(Zo,arguments)),P(this,function(e){e.__transition__[t].ease=n}))},Cs.delay=function(n){var t=this.id;return arguments.length<1?this.node().__transition__[t].delay:P(this,"function"==typeof n?function(e,r,u){e.__transition__[t].delay=+n.call(e,e.__data__,r,u)}:(n=+n,function(e){e.__transition__[t].delay=n}))},Cs.duration=function(n){var t=this.id;return arguments.length<1?this.node().__transition__[t].duration:P(this,"function"==typeof n?function(e,r,u){e.__transition__[t].duration=Math.max(1,n.call(e,e.__data__,r,u))}:(n=Math.max(1,n),function(e){e.__transition__[t].duration=n}))},Cs.each=function(n,t){var e=this.id;if(arguments.length<2){var r=ks,u=Ss;Ss=e,P(this,function(t,r,u){ks=t.__transition__[e],n.call(t,t.__data__,r,u)}),ks=r,Ss=u}else P(this,function(r){var u=r.__transition__[e];(u.event||(u.event=Zo.dispatch("start","end"))).on(n,t)});return this},Cs.transition=function(){for(var n,t,e,r,u=this.id,i=++Ns,o=[],a=0,c=this.length;c>a;a++){o.push(n=[]);for(var t=this[a],s=0,l=t.length;l>s;s++)(e=t[s])&&(r=Object.create(e.__transition__[u]),r.delay+=r.duration,Po(e,s,i,r)),n.push(e)}return qo(o,i)},Zo.svg.axis=function(){function n(n){n.each(function(){var n,s=Zo.select(this),l=this.__chart__||e,f=this.__chart__=e.copy(),h=null==c?f.ticks?f.ticks.apply(f,a):f.domain():c,g=null==t?f.tickFormat?f.tickFormat.apply(f,a):wt:t,p=s.selectAll(".tick").data(h,f),v=p.enter().insert("g",".domain").attr("class","tick").style("opacity",ka),d=Zo.transition(p.exit()).style("opacity",ka).remove(),m=Zo.transition(p.order()).style("opacity",1),y=Ti(f),x=s.selectAll(".domain").data([0]),M=(x.enter().append("path").attr("class","domain"),Zo.transition(x));v.append("line"),v.append("text");var _=v.select("line"),b=m.select("line"),w=p.select("text").text(g),S=v.select("text"),k=m.select("text");switch(r){case"bottom":n=Uo,_.attr("y2",u),S.attr("y",Math.max(u,0)+o),b.attr("x2",0).attr("y2",u),k.attr("x",0).attr("y",Math.max(u,0)+o),w.attr("dy",".71em").style("text-anchor","middle"),M.attr("d","M"+y[0]+","+i+"V0H"+y[1]+"V"+i);break;case"top":n=Uo,_.attr("y2",-u),S.attr("y",-(Math.max(u,0)+o)),b.attr("x2",0).attr("y2",-u),k.attr("x",0).attr("y",-(Math.max(u,0)+o)),w.attr("dy","0em").style("text-anchor","middle"),M.attr("d","M"+y[0]+","+-i+"V0H"+y[1]+"V"+-i);break;case"left":n=jo,_.attr("x2",-u),S.attr("x",-(Math.max(u,0)+o)),b.attr("x2",-u).attr("y2",0),k.attr("x",-(Math.max(u,0)+o)).attr("y",0),w.attr("dy",".32em").style("text-anchor","end"),M.attr("d","M"+-i+","+y[0]+"H0V"+y[1]+"H"+-i);break;case"right":n=jo,_.attr("x2",u),S.attr("x",Math.max(u,0)+o),b.attr("x2",u).attr("y2",0),k.attr("x",Math.max(u,0)+o).attr("y",0),w.attr("dy",".32em").style("text-anchor","start"),M.attr("d","M"+i+","+y[0]+"H0V"+y[1]+"H"+i)}if(f.rangeBand){var E=f,A=E.rangeBand()/2;l=f=function(n){return E(n)+A}}else l.rangeBand?l=f:d.call(n,f);v.call(n,l),m.call(n,f)})}var t,e=Zo.scale.linear(),r=zs,u=6,i=6,o=3,a=[10],c=null;return n.scale=function(t){return arguments.length?(e=t,n):e},n.orient=function(t){return arguments.length?(r=t in Ls?t+"":zs,n):r},n.ticks=function(){return arguments.length?(a=arguments,n):a},n.tickValues=function(t){return arguments.length?(c=t,n):c},n.tickFormat=function(e){return arguments.length?(t=e,n):t},n.tickSize=function(t){var e=arguments.length;return e?(u=+t,i=+arguments[e-1],n):u},n.innerTickSize=function(t){return arguments.length?(u=+t,n):u},n.outerTickSize=function(t){return arguments.length?(i=+t,n):i},n.tickPadding=function(t){return arguments.length?(o=+t,n):o},n.tickSubdivide=function(){return arguments.length&&n},n};var zs="bottom",Ls={top:1,right:1,bottom:1,left:1};Zo.svg.brush=function(){function n(i){i.each(function(){var i=Zo.select(this).style("pointer-events","all").style("-webkit-tap-highlight-color","rgba(0,0,0,0)").on("mousedown.brush",u).on("touchstart.brush",u),o=i.selectAll(".background").data([0]);o.enter().append("rect").attr("class","background").style("visibility","hidden").style("cursor","crosshair"),i.selectAll(".extent").data([0]).enter().append("rect").attr("class","extent").style("cursor","move");var a=i.selectAll(".resize").data(p,wt);a.exit().remove(),a.enter().append("g").attr("class",function(n){return"resize "+n}).style("cursor",function(n){return Ts[n]}).append("rect").attr("x",function(n){return/[ew]$/.test(n)?-3:null}).attr("y",function(n){return/^[ns]/.test(n)?-3:null}).attr("width",6).attr("height",6).style("visibility","hidden"),a.style("display",n.empty()?"none":null);var l,f=Zo.transition(i),h=Zo.transition(o);c&&(l=Ti(c),h.attr("x",l[0]).attr("width",l[1]-l[0]),e(f)),s&&(l=Ti(s),h.attr("y",l[0]).attr("height",l[1]-l[0]),r(f)),t(f)})}function t(n){n.selectAll(".resize").attr("transform",function(n){return"translate("+l[+/e$/.test(n)]+","+f[+/^s/.test(n)]+")"})}function e(n){n.select(".extent").attr("x",l[0]),n.selectAll(".extent,.n>rect,.s>rect").attr("width",l[1]-l[0])}function r(n){n.select(".extent").attr("y",f[0]),n.selectAll(".extent,.e>rect,.w>rect").attr("height",f[1]-f[0])}function u(){function u(){32==Zo.event.keyCode&&(C||(x=null,z[0]-=l[1],z[1]-=f[1],C=2),y())}function p(){32==Zo.event.keyCode&&2==C&&(z[0]+=l[1],z[1]+=f[1],C=0,y())}function v(){var n=Zo.mouse(_),u=!1;M&&(n[0]+=M[0],n[1]+=M[1]),C||(Zo.event.altKey?(x||(x=[(l[0]+l[1])/2,(f[0]+f[1])/2]),z[0]=l[+(n[0]<x[0])],z[1]=f[+(n[1]<x[1])]):x=null),E&&d(n,c,0)&&(e(S),u=!0),A&&d(n,s,1)&&(r(S),u=!0),u&&(t(S),w({type:"brush",mode:C?"move":"resize"}))}function d(n,t,e){var r,u,a=Ti(t),c=a[0],s=a[1],p=z[e],v=e?f:l,d=v[1]-v[0];return C&&(c-=p,s-=d+p),r=(e?g:h)?Math.max(c,Math.min(s,n[e])):n[e],C?u=(r+=p)+d:(x&&(p=Math.max(c,Math.min(s,2*x[e]-r))),r>p?(u=r,r=p):u=p),v[0]!=r||v[1]!=u?(e?o=null:i=null,v[0]=r,v[1]=u,!0):void 0}function m(){v(),S.style("pointer-events","all").selectAll(".resize").style("display",n.empty()?"none":null),Zo.select("body").style("cursor",null),L.on("mousemove.brush",null).on("mouseup.brush",null).on("touchmove.brush",null).on("touchend.brush",null).on("keydown.brush",null).on("keyup.brush",null),N(),w({type:"brushend"})}var x,M,_=this,b=Zo.select(Zo.event.target),w=a.of(_,arguments),S=Zo.select(_),k=b.datum(),E=!/^(n|s)$/.test(k)&&c,A=!/^(e|w)$/.test(k)&&s,C=b.classed("extent"),N=I(),z=Zo.mouse(_),L=Zo.select(Wo).on("keydown.brush",u).on("keyup.brush",p);if(Zo.event.changedTouches?L.on("touchmove.brush",v).on("touchend.brush",m):L.on("mousemove.brush",v).on("mouseup.brush",m),S.interrupt().selectAll("*").interrupt(),C)z[0]=l[0]-z[0],z[1]=f[0]-z[1];else if(k){var T=+/w$/.test(k),q=+/^n/.test(k);M=[l[1-T]-z[0],f[1-q]-z[1]],z[0]=l[T],z[1]=f[q]}else Zo.event.altKey&&(x=z.slice());S.style("pointer-events","none").selectAll(".resize").style("display",null),Zo.select("body").style("cursor",b.style("cursor")),w({type:"brushstart"}),v()}var i,o,a=M(n,"brushstart","brush","brushend"),c=null,s=null,l=[0,0],f=[0,0],h=!0,g=!0,p=qs[0];return n.event=function(n){n.each(function(){var n=a.of(this,arguments),t={x:l,y:f,i:i,j:o},e=this.__chart__||t;this.__chart__=t,Ss?Zo.select(this).transition().each("start.brush",function(){i=e.i,o=e.j,l=e.x,f=e.y,n({type:"brushstart"})}).tween("brush:brush",function(){var e=gu(l,t.x),r=gu(f,t.y);return i=o=null,function(u){l=t.x=e(u),f=t.y=r(u),n({type:"brush",mode:"resize"})}}).each("end.brush",function(){i=t.i,o=t.j,n({type:"brush",mode:"resize"}),n({type:"brushend"})}):(n({type:"brushstart"}),n({type:"brush",mode:"resize"}),n({type:"brushend"}))})},n.x=function(t){return arguments.length?(c=t,p=qs[!c<<1|!s],n):c},n.y=function(t){return arguments.length?(s=t,p=qs[!c<<1|!s],n):s},n.clamp=function(t){return arguments.length?(c&&s?(h=!!t[0],g=!!t[1]):c?h=!!t:s&&(g=!!t),n):c&&s?[h,g]:c?h:s?g:null},n.extent=function(t){var e,r,u,a,h;return arguments.length?(c&&(e=t[0],r=t[1],s&&(e=e[0],r=r[0]),i=[e,r],c.invert&&(e=c(e),r=c(r)),e>r&&(h=e,e=r,r=h),(e!=l[0]||r!=l[1])&&(l=[e,r])),s&&(u=t[0],a=t[1],c&&(u=u[1],a=a[1]),o=[u,a],s.invert&&(u=s(u),a=s(a)),u>a&&(h=u,u=a,a=h),(u!=f[0]||a!=f[1])&&(f=[u,a])),n):(c&&(i?(e=i[0],r=i[1]):(e=l[0],r=l[1],c.invert&&(e=c.invert(e),r=c.invert(r)),e>r&&(h=e,e=r,r=h))),s&&(o?(u=o[0],a=o[1]):(u=f[0],a=f[1],s.invert&&(u=s.invert(u),a=s.invert(a)),u>a&&(h=u,u=a,a=h))),c&&s?[[e,u],[r,a]]:c?[e,r]:s&&[u,a])},n.clear=function(){return n.empty()||(l=[0,0],f=[0,0],i=o=null),n},n.empty=function(){return!!c&&l[0]==l[1]||!!s&&f[0]==f[1]},Zo.rebind(n,a,"on")};var Ts={n:"ns-resize",e:"ew-resize",s:"ns-resize",w:"ew-resize",nw:"nwse-resize",ne:"nesw-resize",se:"nwse-resize",sw:"nesw-resize"},qs=[["n","e","s","w","nw","ne","se","sw"],["e","w"],["n","s"],[]],Rs=Qa.format=ic.timeFormat,Ds=Rs.utc,Ps=Ds("%Y-%m-%dT%H:%M:%S.%LZ");Rs.iso=Date.prototype.toISOString&&+new Date("2000-01-01T00:00:00.000Z")?Ho:Ps,Ho.parse=function(n){var t=new Date(n);return isNaN(t)?null:t},Ho.toString=Ps.toString,Qa.second=Dt(function(n){return new nc(1e3*Math.floor(n/1e3))},function(n,t){n.setTime(n.getTime()+1e3*Math.floor(t))},function(n){return n.getSeconds()}),Qa.seconds=Qa.second.range,Qa.seconds.utc=Qa.second.utc.range,Qa.minute=Dt(function(n){return new nc(6e4*Math.floor(n/6e4))},function(n,t){n.setTime(n.getTime()+6e4*Math.floor(t))},function(n){return n.getMinutes()}),Qa.minutes=Qa.minute.range,Qa.minutes.utc=Qa.minute.utc.range,Qa.hour=Dt(function(n){var t=n.getTimezoneOffset()/60;return new nc(36e5*(Math.floor(n/36e5-t)+t))},function(n,t){n.setTime(n.getTime()+36e5*Math.floor(t))},function(n){return n.getHours()}),Qa.hours=Qa.hour.range,Qa.hours.utc=Qa.hour.utc.range,Qa.month=Dt(function(n){return n=Qa.day(n),n.setDate(1),n},function(n,t){n.setMonth(n.getMonth()+t)},function(n){return n.getMonth()}),Qa.months=Qa.month.range,Qa.months.utc=Qa.month.utc.range;var Us=[1e3,5e3,15e3,3e4,6e4,3e5,9e5,18e5,36e5,108e5,216e5,432e5,864e5,1728e5,6048e5,2592e6,7776e6,31536e6],js=[[Qa.second,1],[Qa.second,5],[Qa.second,15],[Qa.second,30],[Qa.minute,1],[Qa.minute,5],[Qa.minute,15],[Qa.minute,30],[Qa.hour,1],[Qa.hour,3],[Qa.hour,6],[Qa.hour,12],[Qa.day,1],[Qa.day,2],[Qa.week,1],[Qa.month,1],[Qa.month,3],[Qa.year,1]],Hs=Rs.multi([[".%L",function(n){return n.getMilliseconds()}],[":%S",function(n){return n.getSeconds()}],["%I:%M",function(n){return n.getMinutes()}],["%I %p",function(n){return n.getHours()}],["%a %d",function(n){return n.getDay()&&1!=n.getDate()}],["%b %d",function(n){return 1!=n.getDate()}],["%B",function(n){return n.getMonth()}],["%Y",we]]),Fs={range:function(n,t,e){return Zo.range(Math.ceil(n/e)*e,+t,e).map(Oo)},floor:wt,ceil:wt};js.year=Qa.year,Qa.scale=function(){return Fo(Zo.scale.linear(),js,Hs)};var Os=js.map(function(n){return[n[0].utc,n[1]]}),Ys=Ds.multi([[".%L",function(n){return n.getUTCMilliseconds()}],[":%S",function(n){return n.getUTCSeconds()}],["%I:%M",function(n){return n.getUTCMinutes()}],["%I %p",function(n){return n.getUTCHours()}],["%a %d",function(n){return n.getUTCDay()&&1!=n.getUTCDate()}],["%b %d",function(n){return 1!=n.getUTCDate()}],["%B",function(n){return n.getUTCMonth()}],["%Y",we]]);Os.year=Qa.year.utc,Qa.scale.utc=function(){return Fo(Zo.scale.linear(),Os,Ys)},Zo.text=St(function(n){return n.responseText}),Zo.json=function(n,t){return kt(n,"application/json",Yo,t)},Zo.html=function(n,t){return kt(n,"text/html",Io,t)},Zo.xml=St(function(n){return n.responseXML}),"function"==typeof define&&define.amd?define(Zo):"object"==typeof module&&module.exports&&(module.exports=Zo),this.d3=Zo}();
// d3.tip
// Copyright (c) 2013 Justin Palmer
//
// Tooltips for d3.js SVG visualizations

// Public - contructs a new tooltip
//
// Returns a tip
d3.tip = function() {
  var direction = d3_tip_direction,
      offset    = d3_tip_offset,
      html      = d3_tip_html,
      node      = initNode(),
      svg       = null,
      point     = null,
      target    = null

  function tip(vis) {
    svg = getSVGNode(vis)
    point = svg.createSVGPoint()
    document.body.appendChild(node)
  }

  // Public - show the tooltip on the screen
  //
  // Returns a tip
  tip.show = function() {
    var args = Array.prototype.slice.call(arguments)
    if(args[args.length - 1] instanceof SVGElement) target = args.pop()

    var content = html.apply(this, args),
        poffset = offset.apply(this, args),
        dir     = direction.apply(this, args),
        nodel   = d3.select(node), i = 0,
        coords

    nodel.html(content)
      .style({ opacity: 1, 'pointer-events': 'all' })

    while(i--) nodel.classed(directions[i], false)
    coords = direction_callbacks.get(dir).apply(this)
    nodel.classed(dir, true).style({
      top: (coords.top +  poffset[0]) + 'px',
      left: (coords.left + poffset[1]) + 'px'
    })

    return tip
  }

  // Public - hide the tooltip
  //
  // Returns a tip
  tip.hide = function() {
    nodel = d3.select(node)
    nodel.style({ opacity: 0, 'pointer-events': 'none' })
    return tip
  }

  // Public: Proxy attr calls to the d3 tip container.  Sets or gets attribute value.
  //
  // n - name of the attribute
  // v - value of the attribute
  //
  // Returns tip or attribute value
  tip.attr = function(n, v) {
    if (arguments.length < 2 && typeof n === 'string') {
      return d3.select(node).attr(n)
    } else {
      var args =  Array.prototype.slice.call(arguments)
      d3.selection.prototype.attr.apply(d3.select(node), args)
    }

    return tip
  }

  // Public: Proxy style calls to the d3 tip container.  Sets or gets a style value.
  //
  // n - name of the property
  // v - value of the property
  //
  // Returns tip or style property value
  tip.style = function(n, v) {
    if (arguments.length < 2 && typeof n === 'string') {
      return d3.select(node).style(n)
    } else {
      var args =  Array.prototype.slice.call(arguments)
      d3.selection.prototype.style.apply(d3.select(node), args)
    }

    return tip
  }

  // Public: Set or get the direction of the tooltip
  //
  // v - One of n(north), s(south), e(east), or w(west), nw(northwest),
  //     sw(southwest), ne(northeast) or se(southeast)
  //
  // Returns tip or direction
  tip.direction = function(v) {
    if (!arguments.length) return direction
    direction = v == null ? v : d3.functor(v)

    return tip
  }

  // Public: Sets or gets the offset of the tip
  //
  // v - Array of [x, y] offset
  //
  // Returns offset or
  tip.offset = function(v) {
    if (!arguments.length) return offset
    offset = v == null ? v : d3.functor(v)

    return tip
  }

  // Public: sets or gets the html value of the tooltip
  //
  // v - String value of the tip
  //
  // Returns html value or tip
  tip.html = function(v) {
    if (!arguments.length) return html
    html = v == null ? v : d3.functor(v)

    return tip
  }

  function d3_tip_direction() { return 'n' }
  function d3_tip_offset() { return [0, 0] }
  function d3_tip_html() { return ' ' }

  var direction_callbacks = d3.map({
    n:  direction_n,
    s:  direction_s,
    e:  direction_e,
    w:  direction_w,
    nw: direction_nw,
    ne: direction_ne,
    sw: direction_sw,
    se: direction_se
  }),

  directions = direction_callbacks.keys()

  function direction_n() {
    var bbox = getScreenBBox()
    return {
      top:  bbox.n.y - node.offsetHeight,
      left: bbox.n.x - node.offsetWidth / 2
    }
  }

  function direction_s() {
    var bbox = getScreenBBox()
    return {
      top:  bbox.s.y,
      left: bbox.s.x - node.offsetWidth / 2
    }
  }

  function direction_e() {
    var bbox = getScreenBBox()
    return {
      top:  bbox.e.y - node.offsetHeight / 2,
      left: bbox.e.x
    }
  }

  function direction_w() {
    var bbox = getScreenBBox()
    return {
      top:  bbox.w.y - node.offsetHeight / 2,
      left: bbox.w.x - node.offsetWidth
    }
  }

  function direction_nw() {
    var bbox = getScreenBBox()
    return {
      top:  bbox.nw.y - node.offsetHeight,
      left: bbox.nw.x - node.offsetWidth
    }
  }

  function direction_ne() {
    var bbox = getScreenBBox()
    return {
      top:  bbox.ne.y - node.offsetHeight,
      left: bbox.ne.x
    }
  }

  function direction_sw() {
    var bbox = getScreenBBox()
    return {
      top:  bbox.sw.y,
      left: bbox.sw.x - node.offsetWidth
    }
  }

  function direction_se() {
    var bbox = getScreenBBox()
    return {
      top:  bbox.se.y,
      left: bbox.e.x
    }
  }

  function initNode() {
    var node = d3.select(document.createElement('div'))
    node.style({
      position: 'absolute',
      opacity: 0,
      pointerEvents: 'none',
      boxSizing: 'border-box'
    })

    return node.node()
  }

  function getSVGNode(el) {
    el = el.node()
    if(el.tagName.toLowerCase() == 'svg')
      return el

    return el.ownerSVGElement
  }

  // Private - gets the screen coordinates of a shape
  //
  // Given a shape on the screen, will return an SVGPoint for the directions
  // n(north), s(south), e(east), w(west), ne(northeast), se(southeast), nw(northwest),
  // sw(southwest).
  //
  //    +-+-+
  //    |   |
  //    +   +
  //    |   |
  //    +-+-+
  //
  // Returns an Object {n, s, e, w, nw, sw, ne, se}
  function getScreenBBox() {
    var targetel   = target || d3.event.target,
        bbox       = {},
        matrix     = targetel.getScreenCTM(),
        tbbox      = targetel.getBBox(),
        width      = tbbox.width,
        height     = tbbox.height,
        x          = tbbox.x,
        y          = tbbox.y,
        scrollTop  = document.documentElement.scrollTop || document.body.scrollTop,
        scrollLeft = document.documentElement.scrollLeft || document.body.scrollLeft


    point.x = x + scrollLeft
    point.y = y + scrollTop
    bbox.nw = point.matrixTransform(matrix)
    point.x += width
    bbox.ne = point.matrixTransform(matrix)
    point.y += height
    bbox.se = point.matrixTransform(matrix)
    point.x -= width
    bbox.sw = point.matrixTransform(matrix)
    point.y -= height / 2
    bbox.w  = point.matrixTransform(matrix)
    point.x += width
    bbox.e = point.matrixTransform(matrix)
    point.x -= width / 2
    point.y -= height / 2
    bbox.n = point.matrixTransform(matrix)
    point.y += height
    bbox.s = point.matrixTransform(matrix)

    return bbox
  }

  return tip
};

!function(a){"use strict";function b(a){this.owner=a}function c(a,b){if(Object.create)b.prototype=Object.create(a.prototype);else{var c=function(){};c.prototype=a.prototype,b.prototype=new c}return b.prototype.constructor=b,b}function d(a){var b=this.internal=new e(this);b.loadConfig(a),b.init(),function c(a,b,d){Object.keys(a).forEach(function(e){b[e]=a[e].bind(d),Object.keys(a[e]).length>0&&c(a[e],b[e],d)})}(h,this,this)}function e(b){var c=this;c.d3=a.d3?a.d3:"undefined"!=typeof require?require("d3"):void 0,c.api=b,c.config=c.getDefaultConfig(),c.data={},c.cache={},c.axes={}}function f(a){b.call(this,a)}function g(a,b){function c(a,b){a.attr("transform",function(a){return"translate("+Math.ceil(b(a)+u)+", 0)"})}function d(a,b){a.attr("transform",function(a){return"translate(0,"+Math.ceil(b(a))+")"})}function e(a){var b=a[0],c=a[a.length-1];return c>b?[b,c]:[c,b]}function f(a){var b,c,d=[];if(a.ticks)return a.ticks.apply(a,n);for(c=a.domain(),b=Math.ceil(c[0]);b<c[1];b++)d.push(b);return d.length>0&&d[0]>0&&d.unshift(d[0]-(d[1]-d[0])),d}function g(){var a,c=p.copy();return b.isCategory&&(a=p.domain(),c.domain([a[0],a[1]-1])),c}function h(a){var b=m?m(a):a;return"undefined"!=typeof b?b:""}function i(a){if(z)return z;var b={h:11.5,w:5.5};return a.select("text").text(h).each(function(a){var c=this.getBoundingClientRect(),d=h(a),e=c.height,f=d?c.width/d.length:void 0;e&&f&&(b.h=e,b.w=f)}).text(""),z=b,b}function j(c){return b.withoutTransition?c:a.transition(c)}function k(m){m.each(function(){function m(a,c){function d(a,b){f=void 0;for(var h=1;h<b.length;h++)if(" "===b.charAt(h)&&(f=h),e=b.substr(0,h+1),g=U.w*e.length,g>c)return d(a.concat(b.substr(0,f?f:h)),b.slice(f?f+1:h));return a.concat(b)}var e,f,g,i=h(a),j=[];return"[object Array]"===Object.prototype.toString.call(i)?i:((!c||0>=c)&&(c=X?95:b.isCategory?Math.ceil(F(G[1])-F(G[0]))-12:110),d(j,i+""))}function n(a,b){var c=U.h;return 0===b&&(c="left"===q||"right"===q?-((V[a.index]-1)*(U.h/2)-3):".71em"),c}function v(a){var b=p(a)+(o?0:u);return L[0]<b&&b<L[1]?r:0}function w(a){return a?a>0?"start":"end":"middle"}function x(a){return a?"rotate("+a+")":""}function y(a){return a?8*Math.sin(Math.PI*(a/180)):0}function z(a){return a?11.5-2.5*(a/15)*(a>0?1:-1):W}var A,B,C,D=k.g=a.select(this),E=this.__chart__||p,F=this.__chart__=g(),G=t?t:f(F),H=D.selectAll(".tick").data(G,F),I=H.enter().insert("g",".domain").attr("class","tick").style("opacity",1e-6),J=H.exit().remove(),K=j(H).style("opacity",1),L=p.rangeExtent?p.rangeExtent():e(p.range()),M=D.selectAll(".domain").data([0]),N=(M.enter().append("path").attr("class","domain"),j(M));I.append("line"),I.append("text");var O=I.select("line"),P=K.select("line"),Q=I.select("text"),R=K.select("text");b.isCategory?(u=Math.ceil((F(1)-F(0))/2),B=o?0:u,C=o?u:0):u=B=0;var S,T,U=i(D.select(".tick")),V=[],W=Math.max(r,0)+s,X="left"===q||"right"===q;S=H.select("text"),T=S.selectAll("tspan").data(function(a,c){var d=b.tickMultiline?m(a,b.tickWidth):[].concat(h(a));return V[c]=d.length,d.map(function(a){return{index:c,splitted:a}})}),T.enter().append("tspan"),T.exit().remove(),T.text(function(a){return a.splitted});var Y=b.tickTextRotate;switch(q){case"bottom":A=c,O.attr("y2",r),Q.attr("y",W),P.attr("x1",B).attr("x2",B).attr("y2",v),R.attr("x",0).attr("y",z(Y)).style("text-anchor",w(Y)).attr("transform",x(Y)),T.attr("x",0).attr("dy",n).attr("dx",y(Y)),N.attr("d","M"+L[0]+","+l+"V0H"+L[1]+"V"+l);break;case"top":A=c,O.attr("y2",-r),Q.attr("y",-W),P.attr("x2",0).attr("y2",-r),R.attr("x",0).attr("y",-W),S.style("text-anchor","middle"),T.attr("x",0).attr("dy","0em"),N.attr("d","M"+L[0]+","+-l+"V0H"+L[1]+"V"+-l);break;case"left":A=d,O.attr("x2",-r),Q.attr("x",-W),P.attr("x2",-r).attr("y1",C).attr("y2",C),R.attr("x",-W).attr("y",u),S.style("text-anchor","end"),T.attr("x",-W).attr("dy",n),N.attr("d","M"+-l+","+L[0]+"H0V"+L[1]+"H"+-l);break;case"right":A=d,O.attr("x2",r),Q.attr("x",W),P.attr("x2",r).attr("y2",0),R.attr("x",W).attr("y",0),S.style("text-anchor","start"),T.attr("x",W).attr("dy",n),N.attr("d","M"+l+","+L[0]+"H0V"+L[1]+"H"+l)}if(F.rangeBand){var Z=F,$=Z.rangeBand()/2;E=F=function(a){return Z(a)+$}}else E.rangeBand?E=F:J.call(A,F);I.call(A,E),K.call(A,F)})}var l,m,n,o,p=a.scale.linear(),q="bottom",r=6,s=3,t=null,u=0,v=!0;return b=b||{},l=b.withOuterTick?6:0,k.scale=function(a){return arguments.length?(p=a,k):p},k.orient=function(a){return arguments.length?(q=a in{top:1,right:1,bottom:1,left:1}?a+"":"bottom",k):q},k.tickFormat=function(a){return arguments.length?(m=a,k):m},k.tickCentered=function(a){return arguments.length?(o=a,k):o},k.tickOffset=function(){return u},k.tickInterval=function(){var a,c;return b.isCategory?a=2*u:(c=k.g.select("path.domain").node().getTotalLength()-2*l,a=c/k.g.selectAll("line").size()),1/0===a?0:a},k.ticks=function(){return arguments.length?(n=arguments,k):n},k.tickCulling=function(a){return arguments.length?(v=a,k):v},k.tickValues=function(a){if("function"==typeof a)t=function(){return a(p.domain())};else{if(!arguments.length)return t;t=a}return k},k}var h,i,j,k={version:"0.4.10"};k.generate=function(a){return new d(a)},k.chart={fn:d.prototype,internal:{fn:e.prototype,axis:{fn:f.prototype}}},h=k.chart.fn,i=k.chart.internal.fn,j=k.chart.internal.axis.fn,i.init=function(){var a=this,b=a.config;if(a.initParams(),b.data_url)a.convertUrlToData(b.data_url,b.data_mimeType,b.data_keys,a.initWithData);else if(b.data_json)a.initWithData(a.convertJsonToData(b.data_json,b.data_keys));else if(b.data_rows)a.initWithData(a.convertRowsToData(b.data_rows));else{if(!b.data_columns)throw Error("url or json or rows or columns is required.");a.initWithData(a.convertColumnsToData(b.data_columns))}},i.initParams=function(){var a=this,b=a.d3,c=a.config;a.clipId="c3-"+ +new Date+"-clip",a.clipIdForXAxis=a.clipId+"-xaxis",a.clipIdForYAxis=a.clipId+"-yaxis",a.clipIdForGrid=a.clipId+"-grid",a.clipIdForSubchart=a.clipId+"-subchart",a.clipPath=a.getClipPath(a.clipId),a.clipPathForXAxis=a.getClipPath(a.clipIdForXAxis),a.clipPathForYAxis=a.getClipPath(a.clipIdForYAxis),a.clipPathForGrid=a.getClipPath(a.clipIdForGrid),a.clipPathForSubchart=a.getClipPath(a.clipIdForSubchart),a.dragStart=null,a.dragging=!1,a.flowing=!1,a.cancelClick=!1,a.mouseover=!1,a.transiting=!1,a.color=a.generateColor(),a.levelColor=a.generateLevelColor(),a.dataTimeFormat=c.data_xLocaltime?b.time.format:b.time.format.utc,a.axisTimeFormat=c.axis_x_localtime?b.time.format:b.time.format.utc,a.defaultAxisTimeFormat=a.axisTimeFormat.multi([[".%L",function(a){return a.getMilliseconds()}],[":%S",function(a){return a.getSeconds()}],["%I:%M",function(a){return a.getMinutes()}],["%I %p",function(a){return a.getHours()}],["%-m/%-d",function(a){return a.getDay()&&1!==a.getDate()}],["%-m/%-d",function(a){return 1!==a.getDate()}],["%-m/%-d",function(a){return a.getMonth()}],["%Y/%-m/%-d",function(){return!0}]]),a.hiddenTargetIds=[],a.hiddenLegendIds=[],a.focusedTargetIds=[],a.defocusedTargetIds=[],a.xOrient=c.axis_rotated?"left":"bottom",a.yOrient=c.axis_rotated?c.axis_y_inner?"top":"bottom":c.axis_y_inner?"right":"left",a.y2Orient=c.axis_rotated?c.axis_y2_inner?"bottom":"top":c.axis_y2_inner?"left":"right",a.subXOrient=c.axis_rotated?"left":"bottom",a.isLegendRight="right"===c.legend_position,a.isLegendInset="inset"===c.legend_position,a.isLegendTop="top-left"===c.legend_inset_anchor||"top-right"===c.legend_inset_anchor,a.isLegendLeft="top-left"===c.legend_inset_anchor||"bottom-left"===c.legend_inset_anchor,a.legendStep=0,a.legendItemWidth=0,a.legendItemHeight=0,a.currentMaxTickWidths={x:0,y:0,y2:0},a.rotated_padding_left=30,a.rotated_padding_right=c.axis_rotated&&!c.axis_x_show?0:30,a.rotated_padding_top=5,a.withoutFadeIn={},a.intervalForObserveInserted=void 0,a.axes.subx=b.selectAll([])},i.initChartElements=function(){this.initBar&&this.initBar(),this.initLine&&this.initLine(),this.initArc&&this.initArc(),this.initGauge&&this.initGauge(),this.initText&&this.initText()},i.initWithData=function(b){var c,d,e=this,g=e.d3,h=e.config,i=!0;e.axis=new f(e),e.initPie&&e.initPie(),e.initBrush&&e.initBrush(),e.initZoom&&e.initZoom(),e.selectChart=h.bindto?"function"==typeof h.bindto.node?h.bindto:g.select(h.bindto):g.selectAll([]),e.selectChart.empty()&&(e.selectChart=g.select(document.createElement("div")).style("opacity",0),e.observeInserted(e.selectChart),i=!1),e.selectChart.html("").classed("c3",!0),e.data.xs={},e.data.targets=e.convertDataToTargets(b),h.data_filter&&(e.data.targets=e.data.targets.filter(h.data_filter)),h.data_hide&&e.addHiddenTargetIds(h.data_hide===!0?e.mapToIds(e.data.targets):h.data_hide),h.legend_hide&&e.addHiddenLegendIds(h.legend_hide===!0?e.mapToIds(e.data.targets):h.legend_hide),e.hasType("gauge")&&(h.legend_show=!1),e.updateSizes(),e.updateScales(),e.x.domain(g.extent(e.getXDomain(e.data.targets))),e.y.domain(e.getYDomain(e.data.targets,"y")),e.y2.domain(e.getYDomain(e.data.targets,"y2")),e.subX.domain(e.x.domain()),e.subY.domain(e.y.domain()),e.subY2.domain(e.y2.domain()),e.orgXDomain=e.x.domain(),e.brush&&e.brush.scale(e.subX),h.zoom_enabled&&e.zoom.scale(e.x),e.svg=e.selectChart.append("svg").style("overflow","hidden").on("mouseenter",function(){return h.onmouseover.call(e)}).on("mouseleave",function(){return h.onmouseout.call(e)}),c=e.svg.append("defs"),e.clipChart=e.appendClip(c,e.clipId),e.clipXAxis=e.appendClip(c,e.clipIdForXAxis),e.clipYAxis=e.appendClip(c,e.clipIdForYAxis),e.clipGrid=e.appendClip(c,e.clipIdForGrid),e.clipSubchart=e.appendClip(c,e.clipIdForSubchart),e.updateSvgSize(),d=e.main=e.svg.append("g").attr("transform",e.getTranslate("main")),e.initSubchart&&e.initSubchart(),e.initTooltip&&e.initTooltip(),e.initLegend&&e.initLegend(),d.append("text").attr("class",l.text+" "+l.empty).attr("text-anchor","middle").attr("dominant-baseline","middle"),e.initRegion(),e.initGrid(),d.append("g").attr("clip-path",e.clipPath).attr("class",l.chart),h.grid_lines_front&&e.initGridLines(),e.initEventRect(),e.initChartElements(),d.insert("rect",h.zoom_privileged?null:"g."+l.regions).attr("class",l.zoomRect).attr("width",e.width).attr("height",e.height).style("opacity",0).on("dblclick.zoom",null),h.axis_x_extent&&e.brush.extent(e.getDefaultExtent()),e.axis.init(),e.updateTargets(e.data.targets),i&&(e.updateDimension(),e.config.oninit.call(e),e.redraw({withTransition:!1,withTransform:!0,withUpdateXDomain:!0,withUpdateOrgXDomain:!0,withTransitionForAxis:!1})),null==a.onresize&&(a.onresize=e.generateResize()),a.onresize.add&&(a.onresize.add(function(){h.onresize.call(e)}),a.onresize.add(function(){e.api.flush()}),a.onresize.add(function(){h.onresized.call(e)})),e.api.element=e.selectChart.node()},i.smoothLines=function(a,b){var c=this;"grid"===b&&a.each(function(){var a=c.d3.select(this),b=a.attr("x1"),d=a.attr("x2"),e=a.attr("y1"),f=a.attr("y2");a.attr({x1:Math.ceil(b),x2:Math.ceil(d),y1:Math.ceil(e),y2:Math.ceil(f)})})},i.updateSizes=function(){var a=this,b=a.config,c=a.legend?a.getLegendHeight():0,d=a.legend?a.getLegendWidth():0,e=a.isLegendRight||a.isLegendInset?0:c,f=a.hasArcType(),g=b.axis_rotated||f?0:a.getHorizontalAxisHeight("x"),h=b.subchart_show&&!f?b.subchart_size_height+g:0;a.currentWidth=a.getCurrentWidth(),a.currentHeight=a.getCurrentHeight(),a.margin=b.axis_rotated?{top:a.getHorizontalAxisHeight("y2")+a.getCurrentPaddingTop(),right:f?0:a.getCurrentPaddingRight(),bottom:a.getHorizontalAxisHeight("y")+e+a.getCurrentPaddingBottom(),left:h+(f?0:a.getCurrentPaddingLeft())}:{top:4+a.getCurrentPaddingTop(),right:f?0:a.getCurrentPaddingRight(),bottom:g+h+e+a.getCurrentPaddingBottom(),left:f?0:a.getCurrentPaddingLeft()},a.margin2=b.axis_rotated?{top:a.margin.top,right:0/0,bottom:20+e,left:a.rotated_padding_left}:{top:a.currentHeight-h-e,right:0/0,bottom:g+e,left:a.margin.left},a.margin3={top:0,right:0/0,bottom:0,left:0},a.updateSizeForLegend&&a.updateSizeForLegend(c,d),a.width=a.currentWidth-a.margin.left-a.margin.right,a.height=a.currentHeight-a.margin.top-a.margin.bottom,a.width<0&&(a.width=0),a.height<0&&(a.height=0),a.width2=b.axis_rotated?a.margin.left-a.rotated_padding_left-a.rotated_padding_right:a.width,a.height2=b.axis_rotated?a.height:a.currentHeight-a.margin2.top-a.margin2.bottom,a.width2<0&&(a.width2=0),a.height2<0&&(a.height2=0),a.arcWidth=a.width-(a.isLegendRight?d+10:0),a.arcHeight=a.height-(a.isLegendRight?0:10),a.hasType("gauge")&&(a.arcHeight+=a.height-a.getGaugeLabelHeight()),a.updateRadius&&a.updateRadius(),a.isLegendRight&&f&&(a.margin3.left=a.arcWidth/2+1.1*a.radiusExpanded)},i.updateTargets=function(a){var b=this;b.updateTargetsForText(a),b.updateTargetsForBar(a),b.updateTargetsForLine(a),b.hasArcType()&&b.updateTargetsForArc&&b.updateTargetsForArc(a),b.updateTargetsForSubchart&&b.updateTargetsForSubchart(a),b.showTargets()},i.showTargets=function(){var a=this;a.svg.selectAll("."+l.target).filter(function(b){return a.isTargetToShow(b.id)}).transition().duration(a.config.transition_duration).style("opacity",1)},i.redraw=function(a,b){var c,d,e,f,g,h,i,j,k,m,n,o,p,q,r,s,t,u,v,x,y,z,A,B,C,D,E,F,G,H=this,I=H.main,J=H.d3,K=H.config,L=H.getShapeIndices(H.isAreaType),M=H.getShapeIndices(H.isBarType),N=H.getShapeIndices(H.isLineType),O=H.hasArcType(),P=H.filterTargetsToShow(H.data.targets),Q=H.xv.bind(H);if(a=a||{},c=w(a,"withY",!0),d=w(a,"withSubchart",!0),e=w(a,"withTransition",!0),h=w(a,"withTransform",!1),i=w(a,"withUpdateXDomain",!1),j=w(a,"withUpdateOrgXDomain",!1),k=w(a,"withTrimXDomain",!0),p=w(a,"withUpdateXAxis",i),m=w(a,"withLegend",!1),n=w(a,"withEventRect",!0),o=w(a,"withDimension",!0),f=w(a,"withTransitionForExit",e),g=w(a,"withTransitionForAxis",e),v=e?K.transition_duration:0,x=f?v:0,y=g?v:0,b=b||H.axis.generateTransitions(y),m&&K.legend_show?H.updateLegend(H.mapToIds(H.data.targets),a,b):o&&H.updateDimension(!0),H.isCategorized()&&0===P.length&&H.x.domain([0,H.axes.x.selectAll(".tick").size()]),P.length?(H.updateXDomain(P,i,j,k),K.axis_x_tick_values||(B=H.axis.updateXAxisTickValues(P))):(H.xAxis.tickValues([]),H.subXAxis.tickValues([])),K.zoom_rescale&&!a.flow&&(E=H.x.orgDomain()),H.y.domain(H.getYDomain(P,"y",E)),H.y2.domain(H.getYDomain(P,"y2",E)),!K.axis_y_tick_values&&K.axis_y_tick_count&&H.yAxis.tickValues(H.axis.generateTickValues(H.y.domain(),K.axis_y_tick_count)),!K.axis_y2_tick_values&&K.axis_y2_tick_count&&H.y2Axis.tickValues(H.axis.generateTickValues(H.y2.domain(),K.axis_y2_tick_count)),H.axis.redraw(b,O),H.axis.updateLabels(e),(i||p)&&P.length)if(K.axis_x_tick_culling&&B){for(C=1;C<B.length;C++)if(B.length/C<K.axis_x_tick_culling_max){D=C;break}H.svg.selectAll("."+l.axisX+" .tick text").each(function(a){var b=B.indexOf(a);b>=0&&J.select(this).style("display",b%D?"none":"block")})}else H.svg.selectAll("."+l.axisX+" .tick text").style("display","block");q=H.generateDrawArea?H.generateDrawArea(L,!1):void 0,r=H.generateDrawBar?H.generateDrawBar(M):void 0,s=H.generateDrawLine?H.generateDrawLine(N,!1):void 0,t=H.generateXYForText(L,M,N,!0),u=H.generateXYForText(L,M,N,!1),c&&(H.subY.domain(H.getYDomain(P,"y")),H.subY2.domain(H.getYDomain(P,"y2"))),H.tooltip.style("display","none"),H.updateXgridFocus(),I.select("text."+l.text+"."+l.empty).attr("x",H.width/2).attr("y",H.height/2).text(K.data_empty_label_text).transition().style("opacity",P.length?0:1),H.updateGrid(v),H.updateRegion(v),H.updateBar(x),H.updateLine(x),H.updateArea(x),H.updateCircle(),H.hasDataLabel()&&H.updateText(x),H.redrawArc&&H.redrawArc(v,x,h),H.redrawSubchart&&H.redrawSubchart(d,b,v,x,L,M,N),I.selectAll("."+l.selectedCircles).filter(H.isBarType.bind(H)).selectAll("circle").remove(),K.interaction_enabled&&!a.flow&&n&&(H.redrawEventRect(),H.updateZoom&&H.updateZoom()),H.updateCircleY(),F=(H.config.axis_rotated?H.circleY:H.circleX).bind(H),G=(H.config.axis_rotated?H.circleX:H.circleY).bind(H),a.flow&&(A=H.generateFlow({targets:P,flow:a.flow,duration:a.flow.duration,drawBar:r,drawLine:s,drawArea:q,cx:F,cy:G,xv:Q,xForText:t,yForText:u})),(v||A)&&H.isTabVisible()?J.transition().duration(v).each(function(){var b=[];[H.redrawBar(r,!0),H.redrawLine(s,!0),H.redrawArea(q,!0),H.redrawCircle(F,G,!0),H.redrawText(t,u,a.flow,!0),H.redrawRegion(!0),H.redrawGrid(!0)].forEach(function(a){a.forEach(function(a){b.push(a)})}),z=H.generateWait(),b.forEach(function(a){z.add(a)})}).call(z,function(){A&&A(),K.onrendered&&K.onrendered.call(H)}):(H.redrawBar(r),H.redrawLine(s),H.redrawArea(q),H.redrawCircle(F,G),H.redrawText(t,u,a.flow),H.redrawRegion(),H.redrawGrid(),K.onrendered&&K.onrendered.call(H)),H.mapToIds(H.data.targets).forEach(function(a){H.withoutFadeIn[a]=!0})},i.updateAndRedraw=function(a){var b,c=this,d=c.config;a=a||{},a.withTransition=w(a,"withTransition",!0),a.withTransform=w(a,"withTransform",!1),a.withLegend=w(a,"withLegend",!1),a.withUpdateXDomain=!0,a.withUpdateOrgXDomain=!0,a.withTransitionForExit=!1,a.withTransitionForTransform=w(a,"withTransitionForTransform",a.withTransition),c.updateSizes(),a.withLegend&&d.legend_show||(b=c.axis.generateTransitions(a.withTransitionForAxis?d.transition_duration:0),c.updateScales(),c.updateSvgSize(),c.transformAll(a.withTransitionForTransform,b)),c.redraw(a,b)},i.redrawWithoutRescale=function(){this.redraw({withY:!1,withSubchart:!1,withEventRect:!1,withTransitionForAxis:!1})},i.isTimeSeries=function(){return"timeseries"===this.config.axis_x_type},i.isCategorized=function(){return this.config.axis_x_type.indexOf("categor")>=0},i.isCustomX=function(){var a=this,b=a.config;return!a.isTimeSeries()&&(b.data_x||v(b.data_xs))},i.isTimeSeriesY=function(){return"timeseries"===this.config.axis_y_type},i.getTranslate=function(a){var b,c,d=this,e=d.config;return"main"===a?(b=s(d.margin.left),c=s(d.margin.top)):"context"===a?(b=s(d.margin2.left),c=s(d.margin2.top)):"legend"===a?(b=d.margin3.left,c=d.margin3.top):"x"===a?(b=0,c=e.axis_rotated?0:d.height):"y"===a?(b=0,c=e.axis_rotated?d.height:0):"y2"===a?(b=e.axis_rotated?0:d.width,c=e.axis_rotated?1:0):"subx"===a?(b=0,c=e.axis_rotated?0:d.height2):"arc"===a&&(b=d.arcWidth/2,c=d.arcHeight/2),"translate("+b+","+c+")"},i.initialOpacity=function(a){return null!==a.value&&this.withoutFadeIn[a.id]?1:0},i.initialOpacityForCircle=function(a){return null!==a.value&&this.withoutFadeIn[a.id]?this.opacityForCircle(a):0},i.opacityForCircle=function(a){var b=this.config.point_show?1:0;return m(a.value)?this.isScatterType(a)?.5:b:0},i.opacityForText=function(){return this.hasDataLabel()?1:0},i.xx=function(a){return a?this.x(a.x):null},i.xv=function(a){var b=this,c=a.value;return b.isTimeSeries()?c=b.parseDate(a.value):b.isCategorized()&&"string"==typeof a.value&&(c=b.config.axis_x_categories.indexOf(a.value)),Math.ceil(b.x(c))},i.yv=function(a){var b=this,c=a.axis&&"y2"===a.axis?b.y2:b.y;return Math.ceil(c(a.value))},i.subxx=function(a){return a?this.subX(a.x):null},i.transformMain=function(a,b){var c,d,e,f=this;b&&b.axisX?c=b.axisX:(c=f.main.select("."+l.axisX),a&&(c=c.transition())),b&&b.axisY?d=b.axisY:(d=f.main.select("."+l.axisY),a&&(d=d.transition())),b&&b.axisY2?e=b.axisY2:(e=f.main.select("."+l.axisY2),a&&(e=e.transition())),(a?f.main.transition():f.main).attr("transform",f.getTranslate("main")),c.attr("transform",f.getTranslate("x")),d.attr("transform",f.getTranslate("y")),e.attr("transform",f.getTranslate("y2")),f.main.select("."+l.chartArcs).attr("transform",f.getTranslate("arc"))},i.transformAll=function(a,b){var c=this;c.transformMain(a,b),c.config.subchart_show&&c.transformContext(a,b),c.legend&&c.transformLegend(a)},i.updateSvgSize=function(){var a=this,b=a.svg.select(".c3-brush .background");a.svg.attr("width",a.currentWidth).attr("height",a.currentHeight),a.svg.selectAll(["#"+a.clipId,"#"+a.clipIdForGrid]).select("rect").attr("width",a.width).attr("height",a.height),a.svg.select("#"+a.clipIdForXAxis).select("rect").attr("x",a.getXAxisClipX.bind(a)).attr("y",a.getXAxisClipY.bind(a)).attr("width",a.getXAxisClipWidth.bind(a)).attr("height",a.getXAxisClipHeight.bind(a)),a.svg.select("#"+a.clipIdForYAxis).select("rect").attr("x",a.getYAxisClipX.bind(a)).attr("y",a.getYAxisClipY.bind(a)).attr("width",a.getYAxisClipWidth.bind(a)).attr("height",a.getYAxisClipHeight.bind(a)),a.svg.select("#"+a.clipIdForSubchart).select("rect").attr("width",a.width).attr("height",b.size()?b.attr("height"):0),a.svg.select("."+l.zoomRect).attr("width",a.width).attr("height",a.height),a.selectChart.style("max-height",a.currentHeight+"px")},i.updateDimension=function(a){var b=this;a||(b.config.axis_rotated?(b.axes.x.call(b.xAxis),b.axes.subx.call(b.subXAxis)):(b.axes.y.call(b.yAxis),b.axes.y2.call(b.y2Axis))),b.updateSizes(),b.updateScales(),b.updateSvgSize(),b.transformAll(!1)},i.observeInserted=function(b){var c,d=this;return"undefined"==typeof MutationObserver?void a.console.error("MutationObserver not defined."):(c=new MutationObserver(function(e){e.forEach(function(e){"childList"===e.type&&e.previousSibling&&(c.disconnect(),d.intervalForObserveInserted=a.setInterval(function(){b.node().parentNode&&(a.clearInterval(d.intervalForObserveInserted),d.updateDimension(),d.config.oninit.call(d),d.redraw({withTransform:!0,withUpdateXDomain:!0,withUpdateOrgXDomain:!0,withTransition:!1,withTransitionForTransform:!1,withLegend:!0}),b.transition().style("opacity",1))},10))})}),void c.observe(b.node(),{attributes:!0,childList:!0,characterData:!0}))},i.generateResize=function(){function a(){b.forEach(function(a){a()})}var b=[];return a.add=function(a){b.push(a)},a},i.endall=function(a,b){var c=0;a.each(function(){++c}).each("end",function(){--c||b.apply(this,arguments)})},i.generateWait=function(){var a=[],b=function(b,c){var d=setInterval(function(){var b=0;a.forEach(function(a){if(a.empty())return void(b+=1);try{a.transition()}catch(c){b+=1}}),b===a.length&&(clearInterval(d),c&&c())},10)};return b.add=function(b){a.push(b)},b},i.parseDate=function(b){var c,d=this;return b instanceof Date?c=b:"string"==typeof b?c=d.dataTimeFormat(d.config.data_xFormat).parse(b):"number"!=typeof b&&isNaN(b)||(c=new Date(+b)),(!c||isNaN(+c))&&a.console.error("Failed to parse x '"+b+"' to Date object"),c},i.isTabVisible=function(){var a;return"undefined"!=typeof document.hidden?a="hidden":"undefined"!=typeof document.mozHidden?a="mozHidden":"undefined"!=typeof document.msHidden?a="msHidden":"undefined"!=typeof document.webkitHidden&&(a="webkitHidden"),document[a]?!1:!0},i.getDefaultConfig=function(){var a={bindto:"#chart",size_width:void 0,size_height:void 0,padding_left:void 0,padding_right:void 0,padding_top:void 0,padding_bottom:void 0,zoom_enabled:!1,zoom_extent:void 0,zoom_privileged:!1,zoom_rescale:!1,zoom_onzoom:function(){},zoom_onzoomstart:function(){},zoom_onzoomend:function(){},interaction_enabled:!0,onmouseover:function(){},onmouseout:function(){},onresize:function(){},onresized:function(){},oninit:function(){},onrendered:function(){},transition_duration:350,data_x:void 0,data_xs:{},data_xFormat:"%Y-%m-%d",data_xLocaltime:!0,data_xSort:!0,data_idConverter:function(a){return a},data_names:{},data_classes:{},data_groups:[],data_axes:{},data_type:void 0,data_types:{},data_labels:{},data_order:"desc",data_regions:{},data_color:void 0,data_colors:{},data_hide:!1,data_filter:void 0,data_selection_enabled:!1,data_selection_grouped:!1,data_selection_isselectable:function(){return!0},data_selection_multiple:!0,data_selection_draggable:!1,data_onclick:function(){},data_onmouseover:function(){},data_onmouseout:function(){},data_onselected:function(){},data_onunselected:function(){},data_url:void 0,data_json:void 0,data_rows:void 0,data_columns:void 0,data_mimeType:void 0,data_keys:void 0,data_empty_label_text:"",subchart_show:!1,subchart_size_height:60,subchart_onbrush:function(){},color_pattern:[],color_threshold:{},legend_show:!0,legend_hide:!1,legend_position:"bottom",legend_inset_anchor:"top-left",legend_inset_x:10,legend_inset_y:0,legend_inset_step:void 0,legend_item_onclick:void 0,legend_item_onmouseover:void 0,legend_item_onmouseout:void 0,legend_equally:!1,axis_rotated:!1,axis_x_show:!0,axis_x_type:"indexed",axis_x_localtime:!0,axis_x_categories:[],axis_x_tick_centered:!1,axis_x_tick_format:void 0,axis_x_tick_culling:{},axis_x_tick_culling_max:10,axis_x_tick_count:void 0,axis_x_tick_fit:!0,axis_x_tick_values:null,axis_x_tick_rotate:0,axis_x_tick_outer:!0,axis_x_tick_multiline:!0,axis_x_tick_width:null,axis_x_max:void 0,axis_x_min:void 0,axis_x_padding:{},axis_x_height:void 0,axis_x_extent:void 0,axis_x_label:{},axis_y_show:!0,axis_y_type:void 0,axis_y_max:void 0,axis_y_min:void 0,axis_y_inverted:!1,axis_y_center:void 0,axis_y_inner:void 0,axis_y_label:{},axis_y_tick_format:void 0,axis_y_tick_outer:!0,axis_y_tick_values:null,axis_y_tick_count:void 0,axis_y_tick_time_value:void 0,axis_y_tick_time_interval:void 0,axis_y_padding:{},axis_y_default:void 0,axis_y2_show:!1,axis_y2_max:void 0,axis_y2_min:void 0,axis_y2_inverted:!1,axis_y2_center:void 0,axis_y2_inner:void 0,axis_y2_label:{},axis_y2_tick_format:void 0,axis_y2_tick_outer:!0,axis_y2_tick_values:null,axis_y2_tick_count:void 0,axis_y2_padding:{},axis_y2_default:void 0,grid_x_show:!1,grid_x_type:"tick",grid_x_lines:[],grid_y_show:!1,grid_y_lines:[],grid_y_ticks:10,grid_focus_show:!0,grid_lines_front:!0,point_show:!0,point_r:2.5,point_focus_expand_enabled:!0,point_focus_expand_r:void 0,point_select_r:void 0,line_connectNull:!1,line_step_type:"step",bar_width:void 0,bar_width_ratio:.6,bar_width_max:void 0,bar_zerobased:!0,area_zerobased:!0,pie_label_show:!0,pie_label_format:void 0,pie_label_threshold:.05,pie_expand:!0,gauge_label_show:!0,gauge_label_format:void 0,gauge_expand:!0,gauge_min:0,gauge_max:100,gauge_units:void 0,gauge_width:void 0,donut_label_show:!0,donut_label_format:void 0,donut_label_threshold:.05,donut_width:void 0,donut_expand:!0,donut_title:"",regions:[],tooltip_show:!0,tooltip_grouped:!0,tooltip_format_title:void 0,tooltip_format_name:void 0,tooltip_format_value:void 0,tooltip_position:void 0,tooltip_contents:function(a,b,c,d){return this.getTooltipContent?this.getTooltipContent(a,b,c,d):""},tooltip_init_show:!1,tooltip_init_x:0,tooltip_init_position:{top:"0px",left:"50px"}};return Object.keys(this.additionalConfig).forEach(function(b){a[b]=this.additionalConfig[b]},this),a},i.additionalConfig={},i.loadConfig=function(a){function b(){var a=d.shift();return a&&c&&"object"==typeof c&&a in c?(c=c[a],b()):a?void 0:c}var c,d,e,f=this.config;Object.keys(f).forEach(function(g){c=a,d=g.split("_"),e=b(),q(e)&&(f[g]=e)})},i.getScale=function(a,b,c){return(c?this.d3.time.scale():this.d3.scale.linear()).range([a,b])},i.getX=function(a,b,c,d){var e,f=this,g=f.getScale(a,b,f.isTimeSeries()),h=c?g.domain(c):g;f.isCategorized()?(d=d||function(){return 0},g=function(a,b){var c=h(a)+d(a);return b?c:Math.ceil(c)}):g=function(a,b){var c=h(a);return b?c:Math.ceil(c)};for(e in h)g[e]=h[e];return g.orgDomain=function(){return h.domain()},f.isCategorized()&&(g.domain=function(a){return arguments.length?(h.domain(a),g):(a=this.orgDomain(),[a[0],a[1]+1])}),g},i.getY=function(a,b,c){var d=this.getScale(a,b,this.isTimeSeriesY());return c&&d.domain(c),d},i.getYScale=function(a){return"y2"===this.axis.getId(a)?this.y2:this.y},i.getSubYScale=function(a){return"y2"===this.axis.getId(a)?this.subY2:this.subY},i.updateScales=function(){var a=this,b=a.config,c=!a.x;a.xMin=b.axis_rotated?1:0,a.xMax=b.axis_rotated?a.height:a.width,a.yMin=b.axis_rotated?0:a.height,a.yMax=b.axis_rotated?a.width:1,a.subXMin=a.xMin,a.subXMax=a.xMax,a.subYMin=b.axis_rotated?0:a.height2,a.subYMax=b.axis_rotated?a.width2:1,a.x=a.getX(a.xMin,a.xMax,c?void 0:a.x.orgDomain(),function(){return a.xAxis.tickOffset()}),a.y=a.getY(a.yMin,a.yMax,c?b.axis_y_default:a.y.domain()),a.y2=a.getY(a.yMin,a.yMax,c?b.axis_y2_default:a.y2.domain()),a.subX=a.getX(a.xMin,a.xMax,a.orgXDomain,function(b){return b%1?0:a.subXAxis.tickOffset()}),a.subY=a.getY(a.subYMin,a.subYMax,c?b.axis_y_default:a.subY.domain()),a.subY2=a.getY(a.subYMin,a.subYMax,c?b.axis_y2_default:a.subY2.domain()),a.xAxisTickFormat=a.axis.getXAxisTickFormat(),a.xAxisTickValues=a.axis.getXAxisTickValues(),a.yAxisTickValues=a.axis.getYAxisTickValues(),a.y2AxisTickValues=a.axis.getY2AxisTickValues(),a.xAxis=a.axis.getXAxis(a.x,a.xOrient,a.xAxisTickFormat,a.xAxisTickValues,b.axis_x_tick_outer),a.subXAxis=a.axis.getXAxis(a.subX,a.subXOrient,a.xAxisTickFormat,a.xAxisTickValues,b.axis_x_tick_outer),a.yAxis=a.axis.getYAxis(a.y,a.yOrient,b.axis_y_tick_format,a.yAxisTickValues,b.axis_y_tick_outer),a.y2Axis=a.axis.getYAxis(a.y2,a.y2Orient,b.axis_y2_tick_format,a.y2AxisTickValues,b.axis_y2_tick_outer),c||(a.brush&&a.brush.scale(a.subX),b.zoom_enabled&&a.zoom.scale(a.x)),a.updateArc&&a.updateArc()},i.getYDomainMin=function(a){var b,c,d,e,f,g,h=this,i=h.config,j=h.mapToIds(a),k=h.getValuesAsIdKeyed(a);if(i.data_groups.length>0)for(g=h.hasNegativeValueInTargets(a),b=0;b<i.data_groups.length;b++)if(e=i.data_groups[b].filter(function(a){return j.indexOf(a)>=0}),0!==e.length)for(d=e[0],g&&k[d]&&k[d].forEach(function(a,b){k[d][b]=0>a?a:0}),c=1;c<e.length;c++)f=e[c],k[f]&&k[f].forEach(function(a,b){h.axis.getId(f)!==h.axis.getId(d)||!k[d]||g&&+a>0||(k[d][b]+=+a)});return h.d3.min(Object.keys(k).map(function(a){return h.d3.min(k[a])}))},i.getYDomainMax=function(a){var b,c,d,e,f,g,h=this,i=h.config,j=h.mapToIds(a),k=h.getValuesAsIdKeyed(a);if(i.data_groups.length>0)for(g=h.hasPositiveValueInTargets(a),b=0;b<i.data_groups.length;b++)if(e=i.data_groups[b].filter(function(a){return j.indexOf(a)>=0}),0!==e.length)for(d=e[0],g&&k[d]&&k[d].forEach(function(a,b){k[d][b]=a>0?a:0}),c=1;c<e.length;c++)f=e[c],k[f]&&k[f].forEach(function(a,b){h.axis.getId(f)!==h.axis.getId(d)||!k[d]||g&&0>+a||(k[d][b]+=+a)});return h.d3.max(Object.keys(k).map(function(a){return h.d3.max(k[a])}))},i.getYDomain=function(a,b,c){var d,e,f,g,h,i,j,k,l,n,o,p=this,q=p.config,r=a.filter(function(a){return p.axis.getId(a.id)===b}),s=c?p.filterByXDomain(r,c):r,u="y2"===b?q.axis_y2_min:q.axis_y_min,w="y2"===b?q.axis_y2_max:q.axis_y_max,x=p.getYDomainMin(s),y=p.getYDomainMax(s),z="y2"===b?q.axis_y2_center:q.axis_y_center,A=p.hasType("bar",s)&&q.bar_zerobased||p.hasType("area",s)&&q.area_zerobased,B="y2"===b?q.axis_y2_inverted:q.axis_y_inverted,C=p.hasDataLabel()&&q.axis_rotated,D=p.hasDataLabel()&&!q.axis_rotated;return x=m(u)?u:m(w)?w>x?x:w-10:x,y=m(w)?w:m(u)?y>u?y:u+10:y,0===s.length?"y2"===b?p.y2.domain():p.y.domain():(isNaN(x)&&(x=0),isNaN(y)&&(y=x),x===y&&(0>x?y=0:x=0),n=x>=0&&y>=0,o=0>=x&&0>=y,(m(u)&&n||m(w)&&o)&&(A=!1),A&&(n&&(x=0),o&&(y=0)),e=Math.abs(y-x),f=g=h=.1*e,"undefined"!=typeof z&&(i=Math.max(Math.abs(x),Math.abs(y)),y=z+i,x=z-i),C?(j=p.getDataLabelLength(x,y,"width"),k=t(p.y.range()),l=[j[0]/k,j[1]/k],g+=e*(l[1]/(1-l[0]-l[1])),h+=e*(l[0]/(1-l[0]-l[1]))):D&&(j=p.getDataLabelLength(x,y,"height"),g+=p.axis.convertPixelsToAxisPadding(j[1],e),h+=p.axis.convertPixelsToAxisPadding(j[0],e)),"y"===b&&v(q.axis_y_padding)&&(g=p.axis.getPadding(q.axis_y_padding,"top",g,e),h=p.axis.getPadding(q.axis_y_padding,"bottom",h,e)),"y2"===b&&v(q.axis_y2_padding)&&(g=p.axis.getPadding(q.axis_y2_padding,"top",g,e),h=p.axis.getPadding(q.axis_y2_padding,"bottom",h,e)),A&&(n&&(h=x),o&&(g=-y)),d=[x-h,y+g],B?d.reverse():d)},i.getXDomainMin=function(a){var b=this,c=b.config;return q(c.axis_x_min)?b.isTimeSeries()?this.parseDate(c.axis_x_min):c.axis_x_min:b.d3.min(a,function(a){return b.d3.min(a.values,function(a){return a.x})})},i.getXDomainMax=function(a){var b=this,c=b.config;return q(c.axis_x_max)?b.isTimeSeries()?this.parseDate(c.axis_x_max):c.axis_x_max:b.d3.max(a,function(a){return b.d3.max(a.values,function(a){return a.x})})},i.getXDomainPadding=function(a){var b,c,d,e,f=this,g=f.config,h=a[1]-a[0];return f.isCategorized()?c=0:f.hasType("bar")?(b=f.getMaxDataCount(),c=b>1?h/(b-1)/2:.5):c=.01*h,"object"==typeof g.axis_x_padding&&v(g.axis_x_padding)?(d=m(g.axis_x_padding.left)?g.axis_x_padding.left:c,e=m(g.axis_x_padding.right)?g.axis_x_padding.right:c):d=e="number"==typeof g.axis_x_padding?g.axis_x_padding:c,{left:d,right:e}},i.getXDomain=function(a){var b=this,c=[b.getXDomainMin(a),b.getXDomainMax(a)],d=c[0],e=c[1],f=b.getXDomainPadding(c),g=0,h=0;
return d-e!==0||b.isCategorized()||(b.isTimeSeries()?(d=new Date(.5*d.getTime()),e=new Date(1.5*e.getTime())):(d=0===d?1:.5*d,e=0===e?-1:1.5*e)),(d||0===d)&&(g=b.isTimeSeries()?new Date(d.getTime()-f.left):d-f.left),(e||0===e)&&(h=b.isTimeSeries()?new Date(e.getTime()+f.right):e+f.right),[g,h]},i.updateXDomain=function(a,b,c,d,e){var f=this,g=f.config;return c&&(f.x.domain(e?e:f.d3.extent(f.getXDomain(a))),f.orgXDomain=f.x.domain(),g.zoom_enabled&&f.zoom.scale(f.x).updateScaleExtent(),f.subX.domain(f.x.domain()),f.brush&&f.brush.scale(f.subX)),b&&(f.x.domain(e?e:!f.brush||f.brush.empty()?f.orgXDomain:f.brush.extent()),g.zoom_enabled&&f.zoom.scale(f.x).updateScaleExtent()),d&&f.x.domain(f.trimXDomain(f.x.orgDomain())),f.x.domain()},i.trimXDomain=function(a){var b=this;return a[0]<=b.orgXDomain[0]&&(a[1]=+a[1]+(b.orgXDomain[0]-a[0]),a[0]=b.orgXDomain[0]),b.orgXDomain[1]<=a[1]&&(a[0]=+a[0]-(a[1]-b.orgXDomain[1]),a[1]=b.orgXDomain[1]),a},i.isX=function(a){var b=this,c=b.config;return c.data_x&&a===c.data_x||v(c.data_xs)&&x(c.data_xs,a)},i.isNotX=function(a){return!this.isX(a)},i.getXKey=function(a){var b=this,c=b.config;return c.data_x?c.data_x:v(c.data_xs)?c.data_xs[a]:null},i.getXValuesOfXKey=function(a,b){var c,d=this,e=b&&v(b)?d.mapToIds(b):[];return e.forEach(function(b){d.getXKey(b)===a&&(c=d.data.xs[b])}),c},i.getIndexByX=function(a){var b=this,c=b.filterByX(b.data.targets,a);return c.length?c[0].index:null},i.getXValue=function(a,b){var c=this;return a in c.data.xs&&c.data.xs[a]&&m(c.data.xs[a][b])?c.data.xs[a][b]:b},i.getOtherTargetXs=function(){var a=this,b=Object.keys(a.data.xs);return b.length?a.data.xs[b[0]]:null},i.getOtherTargetX=function(a){var b=this.getOtherTargetXs();return b&&a<b.length?b[a]:null},i.addXs=function(a){var b=this;Object.keys(a).forEach(function(c){b.config.data_xs[c]=a[c]})},i.hasMultipleX=function(a){return this.d3.set(Object.keys(a).map(function(b){return a[b]})).size()>1},i.isMultipleX=function(){return v(this.config.data_xs)||!this.config.data_xSort||this.hasType("scatter")},i.addName=function(a){var b,c=this;return a&&(b=c.config.data_names[a.id],a.name=b?b:a.id),a},i.getValueOnIndex=function(a,b){var c=a.filter(function(a){return a.index===b});return c.length?c[0]:null},i.updateTargetX=function(a,b){var c=this;a.forEach(function(a){a.values.forEach(function(d,e){d.x=c.generateTargetX(b[e],a.id,e)}),c.data.xs[a.id]=b})},i.updateTargetXs=function(a,b){var c=this;a.forEach(function(a){b[a.id]&&c.updateTargetX([a],b[a.id])})},i.generateTargetX=function(a,b,c){var d,e=this;return d=e.isTimeSeries()?e.parseDate(a?a:e.getXValue(b,c)):e.isCustomX()&&!e.isCategorized()?m(a)?+a:e.getXValue(b,c):c},i.cloneTarget=function(a){return{id:a.id,id_org:a.id_org,values:a.values.map(function(a){return{x:a.x,value:a.value,id:a.id}})}},i.updateXs=function(){var a=this;a.data.targets.length&&(a.xs=[],a.data.targets[0].values.forEach(function(b){a.xs[b.index]=b.x}))},i.getPrevX=function(a){var b=this.xs[a-1];return"undefined"!=typeof b?b:null},i.getNextX=function(a){var b=this.xs[a+1];return"undefined"!=typeof b?b:null},i.getMaxDataCount=function(){var a=this;return a.d3.max(a.data.targets,function(a){return a.values.length})},i.getMaxDataCountTarget=function(a){var b,c=a.length,d=0;return c>1?a.forEach(function(a){a.values.length>d&&(b=a,d=a.values.length)}):b=c?a[0]:null,b},i.getEdgeX=function(a){var b=this;return a.length?[b.d3.min(a,function(a){return a.values[0].x}),b.d3.max(a,function(a){return a.values[a.values.length-1].x})]:[0,0]},i.mapToIds=function(a){return a.map(function(a){return a.id})},i.mapToTargetIds=function(a){var b=this;return a?o(a)?[a]:a:b.mapToIds(b.data.targets)},i.hasTarget=function(a,b){var c,d=this.mapToIds(a);for(c=0;c<d.length;c++)if(d[c]===b)return!0;return!1},i.isTargetToShow=function(a){return this.hiddenTargetIds.indexOf(a)<0},i.isLegendToShow=function(a){return this.hiddenLegendIds.indexOf(a)<0},i.filterTargetsToShow=function(a){var b=this;return a.filter(function(a){return b.isTargetToShow(a.id)})},i.mapTargetsToUniqueXs=function(a){var b=this,c=b.d3.set(b.d3.merge(a.map(function(a){return a.values.map(function(a){return+a.x})}))).values();return c.map(b.isTimeSeries()?function(a){return new Date(+a)}:function(a){return+a})},i.addHiddenTargetIds=function(a){this.hiddenTargetIds=this.hiddenTargetIds.concat(a)},i.removeHiddenTargetIds=function(a){this.hiddenTargetIds=this.hiddenTargetIds.filter(function(b){return a.indexOf(b)<0})},i.addHiddenLegendIds=function(a){this.hiddenLegendIds=this.hiddenLegendIds.concat(a)},i.removeHiddenLegendIds=function(a){this.hiddenLegendIds=this.hiddenLegendIds.filter(function(b){return a.indexOf(b)<0})},i.getValuesAsIdKeyed=function(a){var b={};return a.forEach(function(a){b[a.id]=[],a.values.forEach(function(c){b[a.id].push(c.value)})}),b},i.checkValueInTargets=function(a,b){var c,d,e,f=Object.keys(a);for(c=0;c<f.length;c++)for(e=a[f[c]].values,d=0;d<e.length;d++)if(b(e[d].value))return!0;return!1},i.hasNegativeValueInTargets=function(a){return this.checkValueInTargets(a,function(a){return 0>a})},i.hasPositiveValueInTargets=function(a){return this.checkValueInTargets(a,function(a){return a>0})},i.isOrderDesc=function(){var a=this.config;return"string"==typeof a.data_order&&"desc"===a.data_order.toLowerCase()},i.isOrderAsc=function(){var a=this.config;return"string"==typeof a.data_order&&"asc"===a.data_order.toLowerCase()},i.orderTargets=function(a){var b=this,c=b.config,d=b.isOrderAsc(),e=b.isOrderDesc();return d||e?a.sort(function(a,b){var c=function(a,b){return a+Math.abs(b.value)},e=a.values.reduce(c,0),f=b.values.reduce(c,0);return d?f-e:e-f}):n(c.data_order)&&a.sort(c.data_order),a},i.filterByX=function(a,b){return this.d3.merge(a.map(function(a){return a.values})).filter(function(a){return a.x-b===0})},i.filterRemoveNull=function(a){return a.filter(function(a){return m(a.value)})},i.filterByXDomain=function(a,b){return a.map(function(a){return{id:a.id,id_org:a.id_org,values:a.values.filter(function(a){return b[0]<=a.x&&a.x<=b[1]})}})},i.hasDataLabel=function(){var a=this.config;return"boolean"==typeof a.data_labels&&a.data_labels?!0:"object"==typeof a.data_labels&&v(a.data_labels)?!0:!1},i.getDataLabelLength=function(a,b,c){var d=this,e=[0,0],f=1.3;return d.selectChart.select("svg").selectAll(".dummy").data([a,b]).enter().append("text").text(function(a){return d.dataLabelFormat(a.id)(a)}).each(function(a,b){e[b]=this.getBoundingClientRect()[c]*f}).remove(),e},i.isNoneArc=function(a){return this.hasTarget(this.data.targets,a.id)},i.isArc=function(a){return"data"in a&&this.hasTarget(this.data.targets,a.data.id)},i.findSameXOfValues=function(a,b){var c,d=a[b].x,e=[];for(c=b-1;c>=0&&d===a[c].x;c--)e.push(a[c]);for(c=b;c<a.length&&d===a[c].x;c++)e.push(a[c]);return e},i.findClosestFromTargets=function(a,b){var c,d=this;return c=a.map(function(a){return d.findClosest(a.values,b)}),d.findClosest(c,b)},i.findClosest=function(a,b){var c,d=this,e=100;return a.filter(function(a){return a&&d.isBarType(a.id)}).forEach(function(a){var b=d.main.select("."+l.bars+d.getTargetSelectorSuffix(a.id)+" ."+l.bar+"-"+a.index).node();!c&&d.isWithinBar(b)&&(c=a)}),a.filter(function(a){return a&&!d.isBarType(a.id)}).forEach(function(a){var f=d.dist(a,b);e>f&&(e=f,c=a)}),c},i.dist=function(a,b){var c=this,d=c.config,e=d.axis_rotated?1:0,f=d.axis_rotated?0:1,g=c.circleY(a,a.index),h=c.x(a.x);return Math.pow(h-b[e],2)+Math.pow(g-b[f],2)},i.convertValuesToStep=function(a){var b,c=[].concat(a);if(!this.isCategorized())return a;for(b=a.length+1;b>0;b--)c[b]=c[b-1];return c[0]={x:c[0].x-1,value:c[0].value,id:c[0].id},c[a.length+1]={x:c[a.length].x+1,value:c[a.length].value,id:c[a.length].id},c},i.updateDataAttributes=function(a,b){var c=this,d=c.config,e=d["data_"+a];return"undefined"==typeof b?e:(Object.keys(b).forEach(function(a){e[a]=b[a]}),c.redraw({withLegend:!0}),e)},i.convertUrlToData=function(a,b,c,d){var e=this,f=b?b:"csv";e.d3.xhr(a,function(a,b){var g;if(!b)throw new Error(a.responseURL+" "+a.status+" ("+a.statusText+")");g="json"===f?e.convertJsonToData(JSON.parse(b.response),c):"tsv"===f?e.convertTsvToData(b.response):e.convertCsvToData(b.response),d.call(e,g)})},i.convertXsvToData=function(a,b){var c,d=b.parseRows(a);return 1===d.length?(c=[{}],d[0].forEach(function(a){c[0][a]=null})):c=b.parse(a),c},i.convertCsvToData=function(a){return this.convertXsvToData(a,this.d3.csv)},i.convertTsvToData=function(a){return this.convertXsvToData(a,this.d3.tsv)},i.convertJsonToData=function(a,b){var c,d,e=this,f=[];return b?(b.x?(c=b.value.concat(b.x),e.config.data_x=b.x):c=b.value,f.push(c),a.forEach(function(a){var b=[];c.forEach(function(c){var d=p(a[c])?null:a[c];b.push(d)}),f.push(b)}),d=e.convertRowsToData(f)):(Object.keys(a).forEach(function(b){f.push([b].concat(a[b]))}),d=e.convertColumnsToData(f)),d},i.convertRowsToData=function(a){var b,c,d=a[0],e={},f=[];for(b=1;b<a.length;b++){for(e={},c=0;c<a[b].length;c++){if(p(a[b][c]))throw new Error("Source data is missing a component at ("+b+","+c+")!");e[d[c]]=a[b][c]}f.push(e)}return f},i.convertColumnsToData=function(a){var b,c,d,e=[];for(b=0;b<a.length;b++)for(d=a[b][0],c=1;c<a[b].length;c++){if(p(e[c-1])&&(e[c-1]={}),p(a[b][c]))throw new Error("Source data is missing a component at ("+b+","+c+")!");e[c-1][d]=a[b][c]}return e},i.convertDataToTargets=function(a,b){var c,d=this,e=d.config,f=d.d3.keys(a[0]).filter(d.isNotX,d),g=d.d3.keys(a[0]).filter(d.isX,d);return f.forEach(function(c){var f=d.getXKey(c);d.isCustomX()||d.isTimeSeries()?g.indexOf(f)>=0?d.data.xs[c]=(b&&d.data.xs[c]?d.data.xs[c]:[]).concat(a.map(function(a){return a[f]}).filter(m).map(function(a,b){return d.generateTargetX(a,c,b)})):e.data_x?d.data.xs[c]=d.getOtherTargetXs():v(e.data_xs)&&(d.data.xs[c]=d.getXValuesOfXKey(f,d.data.targets)):d.data.xs[c]=a.map(function(a,b){return b})}),f.forEach(function(a){if(!d.data.xs[a])throw new Error('x is not defined for id = "'+a+'".')}),c=f.map(function(b,c){var f=e.data_idConverter(b);return{id:f,id_org:b,values:a.map(function(a,g){var h=d.getXKey(b),i=a[h],j=d.generateTargetX(i,b,g);return d.isCustomX()&&d.isCategorized()&&0===c&&i&&(0===g&&(e.axis_x_categories=[]),e.axis_x_categories.push(i)),(p(a[b])||d.data.xs[b].length<=g)&&(j=void 0),{x:j,value:null===a[b]||isNaN(a[b])?null:+a[b],id:f}}).filter(function(a){return q(a.x)})}}),c.forEach(function(a){var b;e.data_xSort&&(a.values=a.values.sort(function(a,b){var c=a.x||0===a.x?a.x:1/0,d=b.x||0===b.x?b.x:1/0;return c-d})),b=0,a.values.forEach(function(a){a.index=b++}),d.data.xs[a.id].sort(function(a,b){return a-b})}),e.data_type&&d.setTargetType(d.mapToIds(c).filter(function(a){return!(a in e.data_types)}),e.data_type),c.forEach(function(a){d.addCache(a.id_org,a)}),c},i.load=function(a,b){var c=this;a&&(b.filter&&(a=a.filter(b.filter)),(b.type||b.types)&&a.forEach(function(a){var d=b.types&&b.types[a.id]?b.types[a.id]:b.type;c.setTargetType(a.id,d)}),c.data.targets.forEach(function(b){for(var c=0;c<a.length;c++)if(b.id===a[c].id){b.values=a[c].values,a.splice(c,1);break}}),c.data.targets=c.data.targets.concat(a)),c.updateTargets(c.data.targets),c.redraw({withUpdateOrgXDomain:!0,withUpdateXDomain:!0,withLegend:!0}),b.done&&b.done()},i.loadFromArgs=function(a){var b=this;a.data?b.load(b.convertDataToTargets(a.data),a):a.url?b.convertUrlToData(a.url,a.mimeType,a.keys,function(c){b.load(b.convertDataToTargets(c),a)}):a.json?b.load(b.convertDataToTargets(b.convertJsonToData(a.json,a.keys)),a):a.rows?b.load(b.convertDataToTargets(b.convertRowsToData(a.rows)),a):a.columns?b.load(b.convertDataToTargets(b.convertColumnsToData(a.columns)),a):b.load(null,a)},i.unload=function(a,b){var c=this;return b||(b=function(){}),a=a.filter(function(a){return c.hasTarget(c.data.targets,a)}),a&&0!==a.length?(c.svg.selectAll(a.map(function(a){return c.selectorTarget(a)})).transition().style("opacity",0).remove().call(c.endall,b),void a.forEach(function(a){c.withoutFadeIn[a]=!1,c.legend&&c.legend.selectAll("."+l.legendItem+c.getTargetSelectorSuffix(a)).remove(),c.data.targets=c.data.targets.filter(function(b){return b.id!==a})})):void b()},i.categoryName=function(a){var b=this.config;return a<b.axis_x_categories.length?b.axis_x_categories[a]:a},i.initEventRect=function(){var a=this;a.main.select("."+l.chart).append("g").attr("class",l.eventRects).style("fill-opacity",0)},i.redrawEventRect=function(){var a,b,c=this,d=c.config,e=c.isMultipleX(),f=c.main.select("."+l.eventRects).style("cursor",d.zoom_enabled?d.axis_rotated?"ns-resize":"ew-resize":null).classed(l.eventRectsMultiple,e).classed(l.eventRectsSingle,!e);f.selectAll("."+l.eventRect).remove(),c.eventRect=f.selectAll("."+l.eventRect),e?(a=c.eventRect.data([0]),c.generateEventRectsForMultipleXs(a.enter()),c.updateEventRect(a)):(b=c.getMaxDataCountTarget(c.data.targets),f.datum(b?b.values:[]),c.eventRect=f.selectAll("."+l.eventRect),a=c.eventRect.data(function(a){return a}),c.generateEventRectsForSingleX(a.enter()),c.updateEventRect(a),a.exit().remove())},i.updateEventRect=function(a){var b,c,d,e,f,g,h=this,i=h.config;a=a||h.eventRect.data(function(a){return a}),h.isMultipleX()?(b=0,c=0,d=h.width,e=h.height):(!h.isCustomX()&&!h.isTimeSeries()||h.isCategorized()?(f=h.getEventRectWidth(),g=function(a){return h.x(a.x)-f/2}):(h.updateXs(),f=function(a){var b=h.getPrevX(a.index),c=h.getNextX(a.index);return null===b&&null===c?i.axis_rotated?h.height:h.width:(null===b&&(b=h.x.domain()[0]),null===c&&(c=h.x.domain()[1]),Math.max(0,(h.x(c)-h.x(b))/2))},g=function(a){var b=h.getPrevX(a.index),c=h.getNextX(a.index),d=h.data.xs[a.id][a.index];return null===b&&null===c?0:(null===b&&(b=h.x.domain()[0]),(h.x(d)+h.x(b))/2)}),b=i.axis_rotated?0:g,c=i.axis_rotated?g:0,d=i.axis_rotated?h.width:f,e=i.axis_rotated?f:h.height),a.attr("class",h.classEvent.bind(h)).attr("x",b).attr("y",c).attr("width",d).attr("height",e)},i.generateEventRectsForSingleX=function(a){var b=this,c=b.d3,d=b.config;a.append("rect").attr("class",b.classEvent.bind(b)).style("cursor",d.data_selection_enabled&&d.data_selection_grouped?"pointer":null).on("mouseover",function(a){var c=a.index;b.dragging||b.flowing||b.hasArcType()||(d.point_focus_expand_enabled&&b.expandCircles(c,null,!0),b.expandBars(c,null,!0),b.main.selectAll("."+l.shape+"-"+c).each(function(a){d.data_onmouseover.call(b.api,a)}))}).on("mouseout",function(a){var c=a.index;b.config&&(b.hasArcType()||(b.hideXGridFocus(),b.hideTooltip(),b.unexpandCircles(),b.unexpandBars(),b.main.selectAll("."+l.shape+"-"+c).each(function(a){d.data_onmouseout.call(b.api,a)})))}).on("mousemove",function(a){var e,f=a.index,g=b.svg.select("."+l.eventRect+"-"+f);b.dragging||b.flowing||b.hasArcType()||(b.isStepType(a)&&"step-after"===b.config.line_step_type&&c.mouse(this)[0]<b.x(b.getXValue(a.id,f))&&(f-=1),e=b.filterTargetsToShow(b.data.targets).map(function(a){return b.addName(b.getValueOnIndex(a.values,f))}),d.tooltip_grouped&&(b.showTooltip(e,this),b.showXGridFocus(e)),(!d.tooltip_grouped||d.data_selection_enabled&&!d.data_selection_grouped)&&b.main.selectAll("."+l.shape+"-"+f).each(function(){c.select(this).classed(l.EXPANDED,!0),d.data_selection_enabled&&g.style("cursor",d.data_selection_grouped?"pointer":null),d.tooltip_grouped||(b.hideXGridFocus(),b.hideTooltip(),d.data_selection_grouped||(b.unexpandCircles(f),b.unexpandBars(f)))}).filter(function(a){return b.isWithinShape(this,a)}).each(function(a){d.data_selection_enabled&&(d.data_selection_grouped||d.data_selection_isselectable(a))&&g.style("cursor","pointer"),d.tooltip_grouped||(b.showTooltip([a],this),b.showXGridFocus([a]),d.point_focus_expand_enabled&&b.expandCircles(f,a.id,!0),b.expandBars(f,a.id,!0))}))}).on("click",function(a){var e=a.index;if(!b.hasArcType()&&b.toggleShape){if(b.cancelClick)return void(b.cancelClick=!1);b.isStepType(a)&&"step-after"===d.line_step_type&&c.mouse(this)[0]<b.x(b.getXValue(a.id,e))&&(e-=1),b.main.selectAll("."+l.shape+"-"+e).each(function(a){(d.data_selection_grouped||b.isWithinShape(this,a))&&(b.toggleShape(this,a,e),b.config.data_onclick.call(b.api,a,this))})}}).call(d.data_selection_draggable&&b.drag?c.behavior.drag().origin(Object).on("drag",function(){b.drag(c.mouse(this))}).on("dragstart",function(){b.dragstart(c.mouse(this))}).on("dragend",function(){b.dragend()}):function(){})},i.generateEventRectsForMultipleXs=function(a){function b(){c.svg.select("."+l.eventRect).style("cursor",null),c.hideXGridFocus(),c.hideTooltip(),c.unexpandCircles(),c.unexpandBars()}var c=this,d=c.d3,e=c.config;a.append("rect").attr("x",0).attr("y",0).attr("width",c.width).attr("height",c.height).attr("class",l.eventRect).on("mouseout",function(){c.config&&(c.hasArcType()||b())}).on("mousemove",function(){var a,f,g,h,i=c.filterTargetsToShow(c.data.targets);if(!c.dragging&&!c.hasArcType(i)){if(a=d.mouse(this),f=c.findClosestFromTargets(i,a),!c.mouseover||f&&f.id===c.mouseover.id||(e.data_onmouseout.call(c.api,c.mouseover),c.mouseover=void 0),!f)return void b();g=c.isScatterType(f)||!e.tooltip_grouped?[f]:c.filterByX(i,f.x),h=g.map(function(a){return c.addName(a)}),c.showTooltip(h,this),e.point_focus_expand_enabled&&c.expandCircles(f.index,f.id,!0),c.expandBars(f.index,f.id,!0),c.showXGridFocus(h),(c.isBarType(f.id)||c.dist(f,a)<100)&&(c.svg.select("."+l.eventRect).style("cursor","pointer"),c.mouseover||(e.data_onmouseover.call(c.api,f),c.mouseover=f))}}).on("click",function(){var a,b,f=c.filterTargetsToShow(c.data.targets);c.hasArcType(f)||(a=d.mouse(this),b=c.findClosestFromTargets(f,a),b&&(c.isBarType(b.id)||c.dist(b,a)<100)&&c.main.selectAll("."+l.shapes+c.getTargetSelectorSuffix(b.id)).selectAll("."+l.shape+"-"+b.index).each(function(){(e.data_selection_grouped||c.isWithinShape(this,b))&&(c.toggleShape(this,b,b.index),c.config.data_onclick.call(c.api,b,this))}))}).call(e.data_selection_draggable&&c.drag?d.behavior.drag().origin(Object).on("drag",function(){c.drag(d.mouse(this))}).on("dragstart",function(){c.dragstart(d.mouse(this))}).on("dragend",function(){c.dragend()}):function(){})},i.dispatchEvent=function(b,c,d){var e=this,f="."+l.eventRect+(e.isMultipleX()?"":"-"+c),g=e.main.select(f).node(),h=g.getBoundingClientRect(),i=h.left+(d?d[0]:0),j=h.top+(d?d[1]:0),k=document.createEvent("MouseEvents");k.initMouseEvent(b,!0,!0,a,0,i,j,i,j,!1,!1,!1,!1,0,null),g.dispatchEvent(k)},i.getCurrentWidth=function(){var a=this,b=a.config;return b.size_width?b.size_width:a.getParentWidth()},i.getCurrentHeight=function(){var a=this,b=a.config,c=b.size_height?b.size_height:a.getParentHeight();return c>0?c:320/(a.hasType("gauge")?2:1)},i.getCurrentPaddingTop=function(){var a=this.config;return m(a.padding_top)?a.padding_top:0},i.getCurrentPaddingBottom=function(){var a=this.config;return m(a.padding_bottom)?a.padding_bottom:0},i.getCurrentPaddingLeft=function(a){var b=this,c=b.config;return m(c.padding_left)?c.padding_left:c.axis_rotated?c.axis_x_show?Math.max(r(b.getAxisWidthByAxisId("x",a)),40):1:!c.axis_y_show||c.axis_y_inner?b.axis.getYAxisLabelPosition().isOuter?30:1:r(b.getAxisWidthByAxisId("y",a))},i.getCurrentPaddingRight=function(){var a=this,b=a.config,c=10,d=a.isLegendRight?a.getLegendWidth()+20:0;return m(b.padding_right)?b.padding_right+1:b.axis_rotated?c+d:!b.axis_y2_show||b.axis_y2_inner?2+d+(a.axis.getY2AxisLabelPosition().isOuter?20:0):r(a.getAxisWidthByAxisId("y2"))+d},i.getParentRectValue=function(a){for(var b,c=this.selectChart.node();c&&"BODY"!==c.tagName;){try{b=c.getBoundingClientRect()[a]}catch(d){"width"===a&&(b=c.offsetWidth)}if(b)break;c=c.parentNode}return b},i.getParentWidth=function(){return this.getParentRectValue("width")},i.getParentHeight=function(){var a=this.selectChart.style("height");return a.indexOf("px")>0?+a.replace("px",""):0},i.getSvgLeft=function(a){var b=this,c=b.config,d=c.axis_rotated||!c.axis_rotated&&!c.axis_y_inner,e=c.axis_rotated?l.axisX:l.axisY,f=b.main.select("."+e).node(),g=f&&d?f.getBoundingClientRect():{right:0},h=b.selectChart.node().getBoundingClientRect(),i=b.hasArcType(),j=g.right-h.left-(i?0:b.getCurrentPaddingLeft(a));return j>0?j:0},i.getAxisWidthByAxisId=function(a,b){var c=this,d=c.axis.getLabelPositionById(a);return c.axis.getMaxTickWidth(a,b)+(d.isInner?20:40)},i.getHorizontalAxisHeight=function(a){var b=this,c=b.config,d=30;return"x"!==a||c.axis_x_show?"x"===a&&c.axis_x_height?c.axis_x_height:"y"!==a||c.axis_y_show?"y2"!==a||c.axis_y2_show?("x"===a&&!c.axis_rotated&&c.axis_x_tick_rotate&&(d=30+b.axis.getMaxTickWidth(a)*Math.cos(Math.PI*(90-c.axis_x_tick_rotate)/180)),d+(b.axis.getLabelPositionById(a).isInner?0:10)+("y2"===a?-10:0)):b.rotated_padding_top:!c.legend_show||b.isLegendRight||b.isLegendInset?1:10:8},i.getEventRectWidth=function(){return Math.max(0,this.xAxis.tickInterval())},i.getShapeIndices=function(a){var b,c,d=this,e=d.config,f={},g=0;return d.filterTargetsToShow(d.data.targets.filter(a,d)).forEach(function(a){for(b=0;b<e.data_groups.length;b++)if(!(e.data_groups[b].indexOf(a.id)<0))for(c=0;c<e.data_groups[b].length;c++)if(e.data_groups[b][c]in f){f[a.id]=f[e.data_groups[b][c]];break}p(f[a.id])&&(f[a.id]=g++)}),f.__max__=g-1,f},i.getShapeX=function(a,b,c,d){var e=this,f=d?e.subX:e.x;return function(d){var e=d.id in c?c[d.id]:0;return d.x||0===d.x?f(d.x)-a*(b/2-e):0}},i.getShapeY=function(a){var b=this;return function(c){var d=a?b.getSubYScale(c.id):b.getYScale(c.id);return d(c.value)}},i.getShapeOffset=function(a,b,c){var d=this,e=d.orderTargets(d.filterTargetsToShow(d.data.targets.filter(a,d))),f=e.map(function(a){return a.id});return function(a,g){var h=c?d.getSubYScale(a.id):d.getYScale(a.id),i=h(0),j=i;return e.forEach(function(c){var e=d.isStepType(a)?d.convertValuesToStep(c.values):c.values;c.id!==a.id&&b[c.id]===b[a.id]&&f.indexOf(c.id)<f.indexOf(a.id)&&e[g].value*a.value>=0&&(j+=h(e[g].value)-i)}),j}},i.isWithinShape=function(a,b){var c,d=this,e=d.d3.select(a);return d.isTargetToShow(b.id)?"circle"===a.nodeName?c=d.isStepType(b)?d.isWithinStep(a,d.getYScale(b.id)(b.value)):d.isWithinCircle(a,1.5*d.pointSelectR(b)):"path"===a.nodeName&&(c=e.classed(l.bar)?d.isWithinBar(a):!0):c=!1,c},i.getInterpolate=function(a){var b=this;return b.isSplineType(a)?"cardinal":b.isStepType(a)?b.config.line_step_type:"linear"},i.initLine=function(){var a=this;a.main.select("."+l.chart).append("g").attr("class",l.chartLines)},i.updateTargetsForLine=function(a){var b,c,d=this,e=d.config,f=d.classChartLine.bind(d),g=d.classLines.bind(d),h=d.classAreas.bind(d),i=d.classCircles.bind(d),j=d.classFocus.bind(d);b=d.main.select("."+l.chartLines).selectAll("."+l.chartLine).data(a).attr("class",function(a){return f(a)+j(a)}),c=b.enter().append("g").attr("class",f).style("opacity",0).style("pointer-events","none"),c.append("g").attr("class",g),c.append("g").attr("class",h),c.append("g").attr("class",function(a){return d.generateClass(l.selectedCircles,a.id)}),c.append("g").attr("class",i).style("cursor",function(a){return e.data_selection_isselectable(a)?"pointer":null}),a.forEach(function(a){d.main.selectAll("."+l.selectedCircles+d.getTargetSelectorSuffix(a.id)).selectAll("."+l.selectedCircle).each(function(b){b.value=a.values[b.index].value})})},i.updateLine=function(a){var b=this;b.mainLine=b.main.selectAll("."+l.lines).selectAll("."+l.line).data(b.lineData.bind(b)),b.mainLine.enter().append("path").attr("class",b.classLine.bind(b)).style("stroke",b.color),b.mainLine.style("opacity",b.initialOpacity.bind(b)).style("shape-rendering",function(a){return b.isStepType(a)?"crispEdges":""}).attr("transform",null),b.mainLine.exit().transition().duration(a).style("opacity",0).remove()},i.redrawLine=function(a,b){return[(b?this.mainLine.transition():this.mainLine).attr("d",a).style("stroke",this.color).style("opacity",1)]},i.generateDrawLine=function(a,b){var c=this,d=c.config,e=c.d3.svg.line(),f=c.generateGetLinePoints(a,b),g=b?c.getSubYScale:c.getYScale,h=function(a){return(b?c.subxx:c.xx).call(c,a)},i=function(a,b){return d.data_groups.length>0?f(a,b)[0][1]:g.call(c,a.id)(a.value)};return e=d.axis_rotated?e.x(i).y(h):e.x(h).y(i),d.line_connectNull||(e=e.defined(function(a){return null!=a.value})),function(a){var f,h=d.line_connectNull?c.filterRemoveNull(a.values):a.values,i=b?c.x:c.subX,j=g.call(c,a.id),k=0,l=0;return c.isLineType(a)?d.data_regions[a.id]?f=c.lineWithRegions(h,i,j,d.data_regions[a.id]):(c.isStepType(a)&&(h=c.convertValuesToStep(h)),f=e.interpolate(c.getInterpolate(a))(h)):(h[0]&&(k=i(h[0].x),l=j(h[0].value)),f=d.axis_rotated?"M "+l+" "+k:"M "+k+" "+l),f?f:"M 0 0"}},i.generateGetLinePoints=function(a,b){var c=this,d=c.config,e=a.__max__+1,f=c.getShapeX(0,e,a,!!b),g=c.getShapeY(!!b),h=c.getShapeOffset(c.isLineType,a,!!b),i=b?c.getSubYScale:c.getYScale;return function(a,b){var e=i.call(c,a.id)(0),j=h(a,b)||e,k=f(a),l=g(a);return d.axis_rotated&&(0<a.value&&e>l||a.value<0&&l>e)&&(l=e),[[k,l-(e-j)],[k,l-(e-j)],[k,l-(e-j)],[k,l-(e-j)]]}},i.lineWithRegions=function(a,b,c,d){function e(a,b){var c;for(c=0;c<b.length;c++)if(b[c].start<a&&a<=b[c].end)return!0;return!1}function f(a){return"M"+a[0][0]+" "+a[0][1]+" "+a[1][0]+" "+a[1][1]}var g,h,i,j,k,l,m,n,o,r,s,t,u=this,v=u.config,w=-1,x="M",y=u.isCategorized()?.5:0,z=[];if(q(d))for(g=0;g<d.length;g++)z[g]={},z[g].start=p(d[g].start)?a[0].x:u.isTimeSeries()?u.parseDate(d[g].start):d[g].start,z[g].end=p(d[g].end)?a[a.length-1].x:u.isTimeSeries()?u.parseDate(d[g].end):d[g].end;for(s=v.axis_rotated?function(a){return c(a.value)}:function(a){return b(a.x)},t=v.axis_rotated?function(a){return b(a.x)}:function(a){return c(a.value)},i=u.isTimeSeries()?function(a,d,e,g){var h,i=a.x.getTime(),j=d.x-a.x,l=new Date(i+j*e),m=new Date(i+j*(e+g));return h=v.axis_rotated?[[c(k(e)),b(l)],[c(k(e+g)),b(m)]]:[[b(l),c(k(e))],[b(m),c(k(e+g))]],f(h)}:function(a,d,e,g){var h;return h=v.axis_rotated?[[c(k(e),!0),b(j(e))],[c(k(e+g),!0),b(j(e+g))]]:[[b(j(e),!0),c(k(e))],[b(j(e+g),!0),c(k(e+g))]],f(h)},g=0;g<a.length;g++){if(p(z)||!e(a[g].x,z))x+=" "+s(a[g])+" "+t(a[g]);else for(j=u.getScale(a[g-1].x+y,a[g].x+y,u.isTimeSeries()),k=u.getScale(a[g-1].value,a[g].value),l=b(a[g].x)-b(a[g-1].x),m=c(a[g].value)-c(a[g-1].value),n=Math.sqrt(Math.pow(l,2)+Math.pow(m,2)),o=2/n,r=2*o,h=o;1>=h;h+=r)x+=i(a[g-1],a[g],h,o);w=a[g].x}return x},i.updateArea=function(a){var b=this,c=b.d3;b.mainArea=b.main.selectAll("."+l.areas).selectAll("."+l.area).data(b.lineData.bind(b)),b.mainArea.enter().append("path").attr("class",b.classArea.bind(b)).style("fill",b.color).style("opacity",function(){return b.orgAreaOpacity=+c.select(this).style("opacity"),0}),b.mainArea.style("opacity",b.orgAreaOpacity),b.mainArea.exit().transition().duration(a).style("opacity",0).remove()},i.redrawArea=function(a,b){return[(b?this.mainArea.transition():this.mainArea).attr("d",a).style("fill",this.color).style("opacity",this.orgAreaOpacity)]},i.generateDrawArea=function(a,b){var c=this,d=c.config,e=c.d3.svg.area(),f=c.generateGetAreaPoints(a,b),g=b?c.getSubYScale:c.getYScale,h=function(a){return(b?c.subxx:c.xx).call(c,a)},i=function(a,b){return d.data_groups.length>0?f(a,b)[0][1]:g.call(c,a.id)(c.getAreaBaseValue(a.id))},j=function(a,b){return d.data_groups.length>0?f(a,b)[1][1]:g.call(c,a.id)(a.value)};return e=d.axis_rotated?e.x0(i).x1(j).y(h):e.x(h).y0(i).y1(j),d.line_connectNull||(e=e.defined(function(a){return null!==a.value})),function(a){var b,f=d.line_connectNull?c.filterRemoveNull(a.values):a.values,g=0,h=0;return c.isAreaType(a)?(c.isStepType(a)&&(f=c.convertValuesToStep(f)),b=e.interpolate(c.getInterpolate(a))(f)):(f[0]&&(g=c.x(f[0].x),h=c.getYScale(a.id)(f[0].value)),b=d.axis_rotated?"M "+h+" "+g:"M "+g+" "+h),b?b:"M 0 0"}},i.getAreaBaseValue=function(){return 0},i.generateGetAreaPoints=function(a,b){var c=this,d=c.config,e=a.__max__+1,f=c.getShapeX(0,e,a,!!b),g=c.getShapeY(!!b),h=c.getShapeOffset(c.isAreaType,a,!!b),i=b?c.getSubYScale:c.getYScale;return function(a,b){var e=i.call(c,a.id)(0),j=h(a,b)||e,k=f(a),l=g(a);return d.axis_rotated&&(0<a.value&&e>l||a.value<0&&l>e)&&(l=e),[[k,j],[k,l-(e-j)],[k,l-(e-j)],[k,j]]}},i.updateCircle=function(){var a=this;a.mainCircle=a.main.selectAll("."+l.circles).selectAll("."+l.circle).data(a.lineOrScatterData.bind(a)),a.mainCircle.enter().append("circle").attr("class",a.classCircle.bind(a)).attr("r",a.pointR.bind(a)).style("fill",a.color),a.mainCircle.style("opacity",a.initialOpacityForCircle.bind(a)),a.mainCircle.exit().remove()},i.redrawCircle=function(a,b,c){var d=this.main.selectAll("."+l.selectedCircle);return[(c?this.mainCircle.transition():this.mainCircle).style("opacity",this.opacityForCircle.bind(this)).style("fill",this.color).attr("cx",a).attr("cy",b),(c?d.transition():d).attr("cx",a).attr("cy",b)]},i.circleX=function(a){return a.x||0===a.x?this.x(a.x):null},i.updateCircleY=function(){var a,b,c=this;c.config.data_groups.length>0?(a=c.getShapeIndices(c.isLineType),b=c.generateGetLinePoints(a),c.circleY=function(a,c){return b(a,c)[0][1]}):c.circleY=function(a){return c.getYScale(a.id)(a.value)}},i.getCircles=function(a,b){var c=this;return(b?c.main.selectAll("."+l.circles+c.getTargetSelectorSuffix(b)):c.main).selectAll("."+l.circle+(m(a)?"-"+a:""))},i.expandCircles=function(a,b,c){var d=this,e=d.pointExpandedR.bind(d);c&&d.unexpandCircles(),d.getCircles(a,b).classed(l.EXPANDED,!0).attr("r",e)},i.unexpandCircles=function(a){var b=this,c=b.pointR.bind(b);b.getCircles(a).filter(function(){return b.d3.select(this).classed(l.EXPANDED)}).classed(l.EXPANDED,!1).attr("r",c)},i.pointR=function(a){var b=this,c=b.config;return b.isStepType(a)?0:n(c.point_r)?c.point_r(a):c.point_r},i.pointExpandedR=function(a){var b=this,c=b.config;return c.point_focus_expand_enabled?c.point_focus_expand_r?c.point_focus_expand_r:1.75*b.pointR(a):b.pointR(a)},i.pointSelectR=function(a){var b=this,c=b.config;return c.point_select_r?c.point_select_r:4*b.pointR(a)},i.isWithinCircle=function(a,b){var c=this.d3,d=c.mouse(a),e=c.select(a),f=+e.attr("cx"),g=+e.attr("cy");return Math.sqrt(Math.pow(f-d[0],2)+Math.pow(g-d[1],2))<b},i.isWithinStep=function(a,b){return Math.abs(b-this.d3.mouse(a)[1])<30},i.initBar=function(){var a=this;a.main.select("."+l.chart).append("g").attr("class",l.chartBars)},i.updateTargetsForBar=function(a){var b,c,d=this,e=d.config,f=d.classChartBar.bind(d),g=d.classBars.bind(d),h=d.classFocus.bind(d);b=d.main.select("."+l.chartBars).selectAll("."+l.chartBar).data(a).attr("class",function(a){return f(a)+h(a)}),c=b.enter().append("g").attr("class",f).style("opacity",0).style("pointer-events","none"),c.append("g").attr("class",g).style("cursor",function(a){return e.data_selection_isselectable(a)?"pointer":null})},i.updateBar=function(a){var b=this,c=b.barData.bind(b),d=b.classBar.bind(b),e=b.initialOpacity.bind(b),f=function(a){return b.color(a.id)};b.mainBar=b.main.selectAll("."+l.bars).selectAll("."+l.bar).data(c),b.mainBar.enter().append("path").attr("class",d).style("stroke",f).style("fill",f),b.mainBar.style("opacity",e),b.mainBar.exit().transition().duration(a).style("opacity",0).remove()},i.redrawBar=function(a,b){return[(b?this.mainBar.transition():this.mainBar).attr("d",a).style("fill",this.color).style("opacity",1)]},i.getBarW=function(a,b){var c=this,d=c.config,e="number"==typeof d.bar_width?d.bar_width:b?a.tickInterval()*d.bar_width_ratio/b:0;return d.bar_width_max&&e>d.bar_width_max?d.bar_width_max:e},i.getBars=function(a,b){var c=this;return(b?c.main.selectAll("."+l.bars+c.getTargetSelectorSuffix(b)):c.main).selectAll("."+l.bar+(m(a)?"-"+a:""))},i.expandBars=function(a,b,c){var d=this;c&&d.unexpandBars(),d.getBars(a,b).classed(l.EXPANDED,!0)},i.unexpandBars=function(a){var b=this;b.getBars(a).classed(l.EXPANDED,!1)},i.generateDrawBar=function(a,b){var c=this,d=c.config,e=c.generateGetBarPoints(a,b);return function(a,b){var c=e(a,b),f=d.axis_rotated?1:0,g=d.axis_rotated?0:1,h="M "+c[0][f]+","+c[0][g]+" L"+c[1][f]+","+c[1][g]+" L"+c[2][f]+","+c[2][g]+" L"+c[3][f]+","+c[3][g]+" z";return h}},i.generateGetBarPoints=function(a,b){var c=this,d=b?c.subXAxis:c.xAxis,e=a.__max__+1,f=c.getBarW(d,e),g=c.getShapeX(f,e,a,!!b),h=c.getShapeY(!!b),i=c.getShapeOffset(c.isBarType,a,!!b),j=b?c.getSubYScale:c.getYScale;
return function(a,b){var d=j.call(c,a.id)(0),e=i(a,b)||d,k=g(a),l=h(a);return c.config.axis_rotated&&(0<a.value&&d>l||a.value<0&&l>d)&&(l=d),[[k,e],[k,l-(d-e)],[k+f,l-(d-e)],[k+f,e]]}},i.isWithinBar=function(a){var b=this.d3.mouse(a),c=a.getBoundingClientRect(),d=a.pathSegList.getItem(0),e=a.pathSegList.getItem(1),f=Math.min(d.x,e.x),g=Math.min(d.y,e.y),h=c.width,i=c.height,j=2,k=f-j,l=f+h+j,m=g+i+j,n=g-j;return k<b[0]&&b[0]<l&&n<b[1]&&b[1]<m},i.initText=function(){var a=this;a.main.select("."+l.chart).append("g").attr("class",l.chartTexts),a.mainText=a.d3.selectAll([])},i.updateTargetsForText=function(a){var b,c,d=this,e=d.classChartText.bind(d),f=d.classTexts.bind(d),g=d.classFocus.bind(d);b=d.main.select("."+l.chartTexts).selectAll("."+l.chartText).data(a).attr("class",function(a){return e(a)+g(a)}),c=b.enter().append("g").attr("class",e).style("opacity",0).style("pointer-events","none"),c.append("g").attr("class",f)},i.updateText=function(a){var b=this,c=b.config,d=b.barOrLineData.bind(b),e=b.classText.bind(b);b.mainText=b.main.selectAll("."+l.texts).selectAll("."+l.text).data(d),b.mainText.enter().append("text").attr("class",e).attr("text-anchor",function(a){return c.axis_rotated?a.value<0?"end":"start":"middle"}).style("stroke","none").style("fill",function(a){return b.color(a)}).style("fill-opacity",0),b.mainText.text(function(a,c,d){return b.dataLabelFormat(a.id)(a.value,a.id,c,d)}),b.mainText.exit().transition().duration(a).style("fill-opacity",0).remove()},i.redrawText=function(a,b,c,d){return[(d?this.mainText.transition():this.mainText).attr("x",a).attr("y",b).style("fill",this.color).style("fill-opacity",c?0:this.opacityForText.bind(this))]},i.getTextRect=function(a,b){var c,d=this.d3.select("body").append("div").classed("c3",!0),e=d.append("svg").style("visibility","hidden").style("position","fixed").style("top",0).style("left",0);return e.selectAll(".dummy").data([a]).enter().append("text").classed(b?b:"",!0).text(a).each(function(){c=this.getBoundingClientRect()}),d.remove(),c},i.generateXYForText=function(a,b,c,d){var e=this,f=e.generateGetAreaPoints(a,!1),g=e.generateGetBarPoints(b,!1),h=e.generateGetLinePoints(c,!1),i=d?e.getXForText:e.getYForText;return function(a,b){var c=e.isAreaType(a)?f:e.isBarType(a)?g:h;return i.call(e,c(a,b),a,this)}},i.getXForText=function(a,b,c){var d,e,f=this,g=c.getBoundingClientRect();return f.config.axis_rotated?(e=f.isBarType(b)?4:6,d=a[2][1]+e*(b.value<0?-1:1)):d=f.hasType("bar")?(a[2][0]+a[0][0])/2:a[0][0],null===b.value&&(d>f.width?d=f.width-g.width:0>d&&(d=4)),d},i.getYForText=function(a,b,c){var d,e=this,f=c.getBoundingClientRect();return e.config.axis_rotated?d=(a[0][0]+a[2][0]+.6*f.height)/2:(d=a[2][1],b.value<0?(d+=f.height,e.isBarType(b)&&e.isSafari()?d-=3:!e.isBarType(b)&&e.isChrome()&&(d+=3)):d+=e.isBarType(b)?-3:-6),null!==b.value||e.config.axis_rotated||(d<f.height?d=f.height:d>this.height&&(d=this.height-4)),d},i.setTargetType=function(a,b){var c=this,d=c.config;c.mapToTargetIds(a).forEach(function(a){c.withoutFadeIn[a]=b===d.data_types[a],d.data_types[a]=b}),a||(d.data_type=b)},i.hasType=function(a,b){var c=this,d=c.config.data_types,e=!1;return b=b||c.data.targets,b&&b.length?b.forEach(function(b){var c=d[b.id];(c&&c.indexOf(a)>=0||!c&&"line"===a)&&(e=!0)}):Object.keys(d).length?Object.keys(d).forEach(function(b){d[b]===a&&(e=!0)}):e=c.config.data_type===a,e},i.hasArcType=function(a){return this.hasType("pie",a)||this.hasType("donut",a)||this.hasType("gauge",a)},i.isLineType=function(a){var b=this.config,c=o(a)?a:a.id;return!b.data_types[c]||["line","spline","area","area-spline","step","area-step"].indexOf(b.data_types[c])>=0},i.isStepType=function(a){var b=o(a)?a:a.id;return["step","area-step"].indexOf(this.config.data_types[b])>=0},i.isSplineType=function(a){var b=o(a)?a:a.id;return["spline","area-spline"].indexOf(this.config.data_types[b])>=0},i.isAreaType=function(a){var b=o(a)?a:a.id;return["area","area-spline","area-step"].indexOf(this.config.data_types[b])>=0},i.isBarType=function(a){var b=o(a)?a:a.id;return"bar"===this.config.data_types[b]},i.isScatterType=function(a){var b=o(a)?a:a.id;return"scatter"===this.config.data_types[b]},i.isPieType=function(a){var b=o(a)?a:a.id;return"pie"===this.config.data_types[b]},i.isGaugeType=function(a){var b=o(a)?a:a.id;return"gauge"===this.config.data_types[b]},i.isDonutType=function(a){var b=o(a)?a:a.id;return"donut"===this.config.data_types[b]},i.isArcType=function(a){return this.isPieType(a)||this.isDonutType(a)||this.isGaugeType(a)},i.lineData=function(a){return this.isLineType(a)?[a]:[]},i.arcData=function(a){return this.isArcType(a.data)?[a]:[]},i.barData=function(a){return this.isBarType(a)?a.values:[]},i.lineOrScatterData=function(a){return this.isLineType(a)||this.isScatterType(a)?a.values:[]},i.barOrLineData=function(a){return this.isBarType(a)||this.isLineType(a)?a.values:[]},i.initGrid=function(){var a=this,b=a.config,c=a.d3;a.grid=a.main.append("g").attr("clip-path",a.clipPathForGrid).attr("class",l.grid),b.grid_x_show&&a.grid.append("g").attr("class",l.xgrids),b.grid_y_show&&a.grid.append("g").attr("class",l.ygrids),b.grid_focus_show&&a.grid.append("g").attr("class",l.xgridFocus).append("line").attr("class",l.xgridFocus),a.xgrid=c.selectAll([]),b.grid_lines_front||a.initGridLines()},i.initGridLines=function(){var a=this,b=a.d3;a.gridLines=a.main.append("g").attr("clip-path",a.clipPathForGrid).attr("class",l.grid+" "+l.gridLines),a.gridLines.append("g").attr("class",l.xgridLines),a.gridLines.append("g").attr("class",l.ygridLines),a.xgridLines=b.selectAll([])},i.updateXGrid=function(a){var b=this,c=b.config,d=b.d3,e=b.generateGridData(c.grid_x_type,b.x),f=b.isCategorized()?b.xAxis.tickOffset():0;b.xgridAttr=c.axis_rotated?{x1:0,x2:b.width,y1:function(a){return b.x(a)-f},y2:function(a){return b.x(a)-f}}:{x1:function(a){return b.x(a)+f},x2:function(a){return b.x(a)+f},y1:0,y2:b.height},b.xgrid=b.main.select("."+l.xgrids).selectAll("."+l.xgrid).data(e),b.xgrid.enter().append("line").attr("class",l.xgrid),a||b.xgrid.attr(b.xgridAttr).style("opacity",function(){return+d.select(this).attr(c.axis_rotated?"y1":"x1")===(c.axis_rotated?b.height:0)?0:1}),b.xgrid.exit().remove()},i.updateYGrid=function(){var a=this,b=a.config,c=a.yAxis.tickValues()||a.y.ticks(b.grid_y_ticks);a.ygrid=a.main.select("."+l.ygrids).selectAll("."+l.ygrid).data(c),a.ygrid.enter().append("line").attr("class",l.ygrid),a.ygrid.attr("x1",b.axis_rotated?a.y:0).attr("x2",b.axis_rotated?a.y:a.width).attr("y1",b.axis_rotated?0:a.y).attr("y2",b.axis_rotated?a.height:a.y),a.ygrid.exit().remove(),a.smoothLines(a.ygrid,"grid")},i.gridTextAnchor=function(a){return a.position?a.position:"end"},i.gridTextDx=function(a){return"start"===a.position?4:"middle"===a.position?0:-4},i.xGridTextX=function(a){return"start"===a.position?-this.height:"middle"===a.position?-this.height/2:0},i.yGridTextX=function(a){return"start"===a.position?0:"middle"===a.position?this.width/2:this.width},i.updateGrid=function(a){var b,c,d,e=this,f=e.main,g=e.config;e.grid.style("visibility",e.hasArcType()?"hidden":"visible"),f.select("line."+l.xgridFocus).style("visibility","hidden"),g.grid_x_show&&e.updateXGrid(),e.xgridLines=f.select("."+l.xgridLines).selectAll("."+l.xgridLine).data(g.grid_x_lines),b=e.xgridLines.enter().append("g").attr("class",function(a){return l.xgridLine+(a["class"]?" "+a["class"]:"")}),b.append("line").style("opacity",0),b.append("text").attr("text-anchor",e.gridTextAnchor).attr("transform",g.axis_rotated?"":"rotate(-90)").attr("dx",e.gridTextDx).attr("dy",-5).style("opacity",0),e.xgridLines.exit().transition().duration(a).style("opacity",0).remove(),g.grid_y_show&&e.updateYGrid(),e.ygridLines=f.select("."+l.ygridLines).selectAll("."+l.ygridLine).data(g.grid_y_lines),c=e.ygridLines.enter().append("g").attr("class",function(a){return l.ygridLine+(a["class"]?" "+a["class"]:"")}),c.append("line").style("opacity",0),c.append("text").attr("text-anchor",e.gridTextAnchor).attr("transform",g.axis_rotated?"rotate(-90)":"").attr("dx",e.gridTextDx).attr("dy",-5).style("opacity",0),d=e.yv.bind(e),e.ygridLines.select("line").transition().duration(a).attr("x1",g.axis_rotated?d:0).attr("x2",g.axis_rotated?d:e.width).attr("y1",g.axis_rotated?0:d).attr("y2",g.axis_rotated?e.height:d).style("opacity",1),e.ygridLines.select("text").transition().duration(a).attr("x",g.axis_rotated?e.xGridTextX.bind(e):e.yGridTextX.bind(e)).attr("y",d).text(function(a){return a.text}).style("opacity",1),e.ygridLines.exit().transition().duration(a).style("opacity",0).remove()},i.redrawGrid=function(a){var b=this,c=b.config,d=b.xv.bind(b),e=b.xgridLines.select("line"),f=b.xgridLines.select("text");return[(a?e.transition():e).attr("x1",c.axis_rotated?0:d).attr("x2",c.axis_rotated?b.width:d).attr("y1",c.axis_rotated?d:0).attr("y2",c.axis_rotated?d:b.height).style("opacity",1),(a?f.transition():f).attr("x",c.axis_rotated?b.yGridTextX.bind(b):b.xGridTextX.bind(b)).attr("y",d).text(function(a){return a.text}).style("opacity",1)]},i.showXGridFocus=function(a){var b=this,c=b.config,d=a.filter(function(a){return a&&m(a.value)}),e=b.main.selectAll("line."+l.xgridFocus),f=b.xx.bind(b);c.tooltip_show&&(b.hasType("scatter")||b.hasArcType()||(e.style("visibility","visible").data([d[0]]).attr(c.axis_rotated?"y1":"x1",f).attr(c.axis_rotated?"y2":"x2",f),b.smoothLines(e,"grid")))},i.hideXGridFocus=function(){this.main.select("line."+l.xgridFocus).style("visibility","hidden")},i.updateXgridFocus=function(){var a=this,b=a.config;a.main.select("line."+l.xgridFocus).attr("x1",b.axis_rotated?0:-10).attr("x2",b.axis_rotated?a.width:-10).attr("y1",b.axis_rotated?-10:0).attr("y2",b.axis_rotated?-10:a.height)},i.generateGridData=function(a,b){var c,d,e,f,g=this,h=[],i=g.main.select("."+l.axisX).selectAll(".tick").size();if("year"===a)for(c=g.getXDomain(),d=c[0].getFullYear(),e=c[1].getFullYear(),f=d;e>=f;f++)h.push(new Date(f+"-01-01 00:00:00"));else h=b.ticks(10),h.length>i&&(h=h.filter(function(a){return(""+a).indexOf(".")<0}));return h},i.getGridFilterToRemove=function(a){return a?function(b){var c=!1;return[].concat(a).forEach(function(a){("value"in a&&b.value===a.value||"class"in a&&b["class"]===a["class"])&&(c=!0)}),c}:function(){return!0}},i.removeGridLines=function(a,b){var c=this,d=c.config,e=c.getGridFilterToRemove(a),f=function(a){return!e(a)},g=b?l.xgridLines:l.ygridLines,h=b?l.xgridLine:l.ygridLine;c.main.select("."+g).selectAll("."+h).filter(e).transition().duration(d.transition_duration).style("opacity",0).remove(),b?d.grid_x_lines=d.grid_x_lines.filter(f):d.grid_y_lines=d.grid_y_lines.filter(f)},i.initTooltip=function(){var a,b=this,c=b.config;if(b.tooltip=b.selectChart.style("position","relative").append("div").attr("class",l.tooltipContainer).style("position","absolute").style("pointer-events","none").style("display","none"),c.tooltip_init_show){if(b.isTimeSeries()&&o(c.tooltip_init_x)){for(c.tooltip_init_x=b.parseDate(c.tooltip_init_x),a=0;a<b.data.targets[0].values.length&&b.data.targets[0].values[a].x-c.tooltip_init_x!==0;a++);c.tooltip_init_x=a}b.tooltip.html(c.tooltip_contents.call(b,b.data.targets.map(function(a){return b.addName(a.values[c.tooltip_init_x])}),b.axis.getXAxisTickFormat(),b.getYFormat(b.hasArcType()),b.color)),b.tooltip.style("top",c.tooltip_init_position.top).style("left",c.tooltip_init_position.left).style("display","block")}},i.getTooltipContent=function(a,b,c,d){var e,f,g,h,i,j,k=this,m=k.config,n=m.tooltip_format_title||b,o=m.tooltip_format_name||function(a){return a},p=m.tooltip_format_value||c;for(f=0;f<a.length;f++)a[f]&&(a[f].value||0===a[f].value)&&(e||(g=n?n(a[f].x):a[f].x,e="<table class='"+l.tooltip+"'>"+(g||0===g?"<tr><th colspan='2'>"+g+"</th></tr>":"")),h=p(a[f].value,a[f].ratio,a[f].id,a[f].index),void 0!==h&&(i=o(a[f].name,a[f].ratio,a[f].id,a[f].index),j=k.levelColor?k.levelColor(a[f].value):d(a[f].id),e+="<tr class='"+l.tooltipName+"-"+a[f].id+"'>",e+="<td class='name'><span style='background-color:"+j+"'></span>"+i+"</td>",e+="<td class='value'>"+h+"</td>",e+="</tr>"));return e+"</table>"},i.tooltipPosition=function(a,b,c,d){var e,f,g,h,i,j=this,k=j.config,l=j.d3,m=j.hasArcType(),n=l.mouse(d);return m?(f=(j.width-(j.isLegendRight?j.getLegendWidth():0))/2+n[0],h=j.height/2+n[1]+20):(e=j.getSvgLeft(!0),k.axis_rotated?(f=e+n[0]+100,g=f+b,i=j.currentWidth-j.getCurrentPaddingRight(),h=j.x(a[0].x)+20):(f=e+j.getCurrentPaddingLeft(!0)+j.x(a[0].x)+20,g=f+b,i=e+j.currentWidth-j.getCurrentPaddingRight(),h=n[1]+15),g>i&&(f-=g-i+20),h+c>j.currentHeight&&(h-=c+30)),0>h&&(h=0),{top:h,left:f}},i.showTooltip=function(a,b){var c,d,e,f=this,g=f.config,h=f.hasArcType(),j=a.filter(function(a){return a&&m(a.value)}),k=g.tooltip_position||i.tooltipPosition;0!==j.length&&g.tooltip_show&&(f.tooltip.html(g.tooltip_contents.call(f,a,f.axis.getXAxisTickFormat(),f.getYFormat(h),f.color)).style("display","block"),c=f.tooltip.property("offsetWidth"),d=f.tooltip.property("offsetHeight"),e=k.call(this,j,c,d,b),f.tooltip.style("top",e.top+"px").style("left",e.left+"px"))},i.hideTooltip=function(){this.tooltip.style("display","none")},i.initLegend=function(){var a=this;return a.legendItemTextBox={},a.legendHasRendered=!1,a.legend=a.svg.append("g").attr("transform",a.getTranslate("legend")),a.config.legend_show?void a.updateLegendWithDefaults():(a.legend.style("visibility","hidden"),void(a.hiddenLegendIds=a.mapToIds(a.data.targets)))},i.updateLegendWithDefaults=function(){var a=this;a.updateLegend(a.mapToIds(a.data.targets),{withTransform:!1,withTransitionForTransform:!1,withTransition:!1})},i.updateSizeForLegend=function(a,b){var c=this,d=c.config,e={top:c.isLegendTop?c.getCurrentPaddingTop()+d.legend_inset_y+5.5:c.currentHeight-a-c.getCurrentPaddingBottom()-d.legend_inset_y,left:c.isLegendLeft?c.getCurrentPaddingLeft()+d.legend_inset_x+.5:c.currentWidth-b-c.getCurrentPaddingRight()-d.legend_inset_x+.5};c.margin3={top:c.isLegendRight?0:c.isLegendInset?e.top:c.currentHeight-a,right:0/0,bottom:0,left:c.isLegendRight?c.currentWidth-b:c.isLegendInset?e.left:0}},i.transformLegend=function(a){var b=this;(a?b.legend.transition():b.legend).attr("transform",b.getTranslate("legend"))},i.updateLegendStep=function(a){this.legendStep=a},i.updateLegendItemWidth=function(a){this.legendItemWidth=a},i.updateLegendItemHeight=function(a){this.legendItemHeight=a},i.getLegendWidth=function(){var a=this;return a.config.legend_show?a.isLegendRight||a.isLegendInset?a.legendItemWidth*(a.legendStep+1):a.currentWidth:0},i.getLegendHeight=function(){var a=this,b=0;return a.config.legend_show&&(b=a.isLegendRight?a.currentHeight:Math.max(20,a.legendItemHeight)*(a.legendStep+1)),b},i.opacityForLegend=function(a){return a.classed(l.legendItemHidden)?null:1},i.opacityForUnfocusedLegend=function(a){return a.classed(l.legendItemHidden)?null:.3},i.toggleFocusLegend=function(a,b){var c=this;a=c.mapToTargetIds(a),c.legend.selectAll("."+l.legendItem).filter(function(b){return a.indexOf(b)>=0}).classed(l.legendItemFocused,b).transition().duration(100).style("opacity",function(){var a=b?c.opacityForLegend:c.opacityForUnfocusedLegend;return a.call(c,c.d3.select(this))})},i.revertLegend=function(){var a=this,b=a.d3;a.legend.selectAll("."+l.legendItem).classed(l.legendItemFocused,!1).transition().duration(100).style("opacity",function(){return a.opacityForLegend(b.select(this))})},i.showLegend=function(a){var b=this,c=b.config;c.legend_show||(c.legend_show=!0,b.legend.style("visibility","visible"),b.legendHasRendered||b.updateLegendWithDefaults()),b.removeHiddenLegendIds(a),b.legend.selectAll(b.selectorLegends(a)).style("visibility","visible").transition().style("opacity",function(){return b.opacityForLegend(b.d3.select(this))})},i.hideLegend=function(a){var b=this,c=b.config;c.legend_show&&u(a)&&(c.legend_show=!1,b.legend.style("visibility","hidden")),b.addHiddenLegendIds(a),b.legend.selectAll(b.selectorLegends(a)).style("opacity",0).style("visibility","hidden")},i.clearLegendItemTextBoxCache=function(){this.legendItemTextBox={}},i.updateLegend=function(a,b,c){function d(a,b){return u.legendItemTextBox[b]||(u.legendItemTextBox[b]=u.getTextRect(a.textContent,l.legendItem)),u.legendItemTextBox[b]}function e(b,c,e){function f(a,b){b||(g=(o-D-n)/2,B>g&&(g=(o-n)/2,D=0,J++)),I[a]=J,H[J]=u.isLegendInset?10:g,E[a]=D,D+=n}var g,h,i=0===e,j=e===a.length-1,k=d(b,c),l=k.width+C+(!j||u.isLegendRight||u.isLegendInset?y:0),m=k.height+x,n=u.isLegendRight||u.isLegendInset?m:l,o=u.isLegendRight||u.isLegendInset?u.getLegendHeight():u.getLegendWidth();return i&&(D=0,J=0,z=0,A=0),v.legend_show&&!u.isLegendToShow(c)?void(F[c]=G[c]=I[c]=E[c]=0):(F[c]=l,G[c]=m,(!z||l>=z)&&(z=l),(!A||m>=A)&&(A=m),h=u.isLegendRight||u.isLegendInset?A:z,void(v.legend_equally?(Object.keys(F).forEach(function(a){F[a]=z}),Object.keys(G).forEach(function(a){G[a]=A}),g=(o-h*a.length)/2,B>g?(D=0,J=0,a.forEach(function(a){f(a)})):f(c,!0)):f(c)))}var f,g,h,i,j,k,m,n,o,p,r,s,t,u=this,v=u.config,x=4,y=10,z=0,A=0,B=10,C=15,D=0,E={},F={},G={},H=[0],I={},J=0;b=b||{},n=w(b,"withTransition",!0),o=w(b,"withTransitionForTransform",!0),u.isLegendInset&&(J=v.legend_inset_step?v.legend_inset_step:a.length,u.updateLegendStep(J)),u.isLegendRight?(f=function(a){return z*I[a]},i=function(a){return H[I[a]]+E[a]}):u.isLegendInset?(f=function(a){return z*I[a]+10},i=function(a){return H[I[a]]+E[a]}):(f=function(a){return H[I[a]]+E[a]},i=function(a){return A*I[a]}),g=function(a,b){return f(a,b)+14},j=function(a,b){return i(a,b)+9},h=function(a,b){return f(a,b)},k=function(a,b){return i(a,b)-5},m=u.legend.selectAll("."+l.legendItem).data(a).enter().append("g").attr("class",function(a){return u.generateClass(l.legendItem,a)}).style("visibility",function(a){return u.isLegendToShow(a)?"visible":"hidden"}).style("cursor","pointer").on("click",function(a){v.legend_item_onclick?v.legend_item_onclick.call(u,a):u.d3.event.altKey?(u.api.hide(),u.api.show(a)):(u.api.toggle(a),u.isTargetToShow(a)?u.api.focus(a):u.api.revert())}).on("mouseover",function(a){u.d3.select(this).classed(l.legendItemFocused,!0),!u.transiting&&u.isTargetToShow(a)&&u.api.focus(a),v.legend_item_onmouseover&&v.legend_item_onmouseover.call(u,a)}).on("mouseout",function(a){u.d3.select(this).classed(l.legendItemFocused,!1),u.api.revert(),v.legend_item_onmouseout&&v.legend_item_onmouseout.call(u,a)}),m.append("text").text(function(a){return q(v.data_names[a])?v.data_names[a]:a}).each(function(a,b){e(this,a,b)}).style("pointer-events","none").attr("x",u.isLegendRight||u.isLegendInset?g:-200).attr("y",u.isLegendRight||u.isLegendInset?-200:j),m.append("rect").attr("class",l.legendItemEvent).style("fill-opacity",0).attr("x",u.isLegendRight||u.isLegendInset?h:-200).attr("y",u.isLegendRight||u.isLegendInset?-200:k),m.append("rect").attr("class",l.legendItemTile).style("pointer-events","none").style("fill",u.color).attr("x",u.isLegendRight||u.isLegendInset?g:-200).attr("y",u.isLegendRight||u.isLegendInset?-200:i).attr("width",10).attr("height",10),t=u.legend.select("."+l.legendBackground+" rect"),u.isLegendInset&&z>0&&0===t.size()&&(t=u.legend.insert("g","."+l.legendItem).attr("class",l.legendBackground).append("rect")),p=u.legend.selectAll("text").data(a).text(function(a){return q(v.data_names[a])?v.data_names[a]:a}).each(function(a,b){e(this,a,b)}),(n?p.transition():p).attr("x",g).attr("y",j),r=u.legend.selectAll("rect."+l.legendItemEvent).data(a),(n?r.transition():r).attr("width",function(a){return F[a]}).attr("height",function(a){return G[a]}).attr("x",h).attr("y",k),s=u.legend.selectAll("rect."+l.legendItemTile).data(a),(n?s.transition():s).style("fill",u.color).attr("x",f).attr("y",i),t&&(n?t.transition():t).attr("height",u.getLegendHeight()-12).attr("width",z*(J+1)+10),u.legend.selectAll("."+l.legendItem).classed(l.legendItemHidden,function(a){return!u.isTargetToShow(a)}),u.updateLegendItemWidth(z),u.updateLegendItemHeight(A),u.updateLegendStep(J),u.updateSizes(),u.updateScales(),u.updateSvgSize(),u.transformAll(o,c),u.legendHasRendered=!0},c(b,f),f.prototype.init=function(){var a=this.owner,b=a.config,c=a.main;a.axes.x=c.append("g").attr("class",l.axis+" "+l.axisX).attr("clip-path",a.clipPathForXAxis).attr("transform",a.getTranslate("x")).style("visibility",b.axis_x_show?"visible":"hidden"),a.axes.x.append("text").attr("class",l.axisXLabel).attr("transform",b.axis_rotated?"rotate(-90)":"").style("text-anchor",this.textAnchorForXAxisLabel.bind(this)),a.axes.y=c.append("g").attr("class",l.axis+" "+l.axisY).attr("clip-path",b.axis_y_inner?"":a.clipPathForYAxis).attr("transform",a.getTranslate("y")).style("visibility",b.axis_y_show?"visible":"hidden"),a.axes.y.append("text").attr("class",l.axisYLabel).attr("transform",b.axis_rotated?"":"rotate(-90)").style("text-anchor",this.textAnchorForYAxisLabel.bind(this)),a.axes.y2=c.append("g").attr("class",l.axis+" "+l.axisY2).attr("transform",a.getTranslate("y2")).style("visibility",b.axis_y2_show?"visible":"hidden"),a.axes.y2.append("text").attr("class",l.axisY2Label).attr("transform",b.axis_rotated?"":"rotate(-90)").style("text-anchor",this.textAnchorForY2AxisLabel.bind(this))},f.prototype.getXAxis=function(a,b,c,d,e,f,h){var i=this.owner,j=i.config,k={isCategory:i.isCategorized(),withOuterTick:e,tickMultiline:j.axis_x_tick_multiline,tickWidth:j.axis_x_tick_width,tickTextRotate:h?0:j.axis_x_tick_rotate,withoutTransition:f},l=g(i.d3,k).scale(a).orient(b);return i.isTimeSeries()&&d&&(d=d.map(function(a){return i.parseDate(a)})),l.tickFormat(c).tickValues(d),i.isCategorized()&&(l.tickCentered(j.axis_x_tick_centered),u(j.axis_x_tick_culling)&&(j.axis_x_tick_culling=!1)),l},f.prototype.updateXAxisTickValues=function(a,b){var c,d=this.owner,e=d.config;return(e.axis_x_tick_fit||e.axis_x_tick_count)&&(c=this.generateTickValues(d.mapTargetsToUniqueXs(a),e.axis_x_tick_count,d.isTimeSeries())),b?b.tickValues(c):(d.xAxis.tickValues(c),d.subXAxis.tickValues(c)),c},f.prototype.getYAxis=function(a,b,c,d,e,f){var h={withOuterTick:e,withoutTransition:f},i=this.owner,j=i.d3,k=i.config,l=g(j,h).scale(a).orient(b).tickFormat(c);return i.isTimeSeriesY()?l.ticks(j.time[k.axis_y_tick_time_value],k.axis_y_tick_time_interval):l.tickValues(d),l},f.prototype.getId=function(a){var b=this.owner.config;return a in b.data_axes?b.data_axes[a]:"y"},f.prototype.getXAxisTickFormat=function(){var a=this.owner,b=a.config,c=a.isTimeSeries()?a.defaultAxisTimeFormat:a.isCategorized()?a.categoryName:function(a){return 0>a?a.toFixed(0):a};return b.axis_x_tick_format&&(n(b.axis_x_tick_format)?c=b.axis_x_tick_format:a.isTimeSeries()&&(c=function(c){return c?a.axisTimeFormat(b.axis_x_tick_format)(c):""})),n(c)?function(b){return c.call(a,b)}:c},f.prototype.getTickValues=function(a,b){return a?a:b?b.tickValues():void 0},f.prototype.getXAxisTickValues=function(){return this.getTickValues(this.owner.config.axis_x_tick_values,this.owner.xAxis)},f.prototype.getYAxisTickValues=function(){return this.getTickValues(this.owner.config.axis_y_tick_values,this.owner.yAxis)},f.prototype.getY2AxisTickValues=function(){return this.getTickValues(this.owner.config.axis_y2_tick_values,this.owner.y2Axis)},f.prototype.getLabelOptionByAxisId=function(a){var b,c=this.owner,d=c.config;return"y"===a?b=d.axis_y_label:"y2"===a?b=d.axis_y2_label:"x"===a&&(b=d.axis_x_label),b},f.prototype.getLabelText=function(a){var b=this.getLabelOptionByAxisId(a);return o(b)?b:b?b.text:null},f.prototype.setLabelText=function(a,b){var c=this.owner,d=c.config,e=this.getLabelOptionByAxisId(a);o(e)?"y"===a?d.axis_y_label=b:"y2"===a?d.axis_y2_label=b:"x"===a&&(d.axis_x_label=b):e&&(e.text=b)},f.prototype.getLabelPosition=function(a,b){var c=this.getLabelOptionByAxisId(a),d=c&&"object"==typeof c&&c.position?c.position:b;return{isInner:d.indexOf("inner")>=0,isOuter:d.indexOf("outer")>=0,isLeft:d.indexOf("left")>=0,isCenter:d.indexOf("center")>=0,isRight:d.indexOf("right")>=0,isTop:d.indexOf("top")>=0,isMiddle:d.indexOf("middle")>=0,isBottom:d.indexOf("bottom")>=0}},f.prototype.getXAxisLabelPosition=function(){return this.getLabelPosition("x",this.owner.config.axis_rotated?"inner-top":"inner-right")},f.prototype.getYAxisLabelPosition=function(){return this.getLabelPosition("y",this.owner.config.axis_rotated?"inner-right":"inner-top")},f.prototype.getY2AxisLabelPosition=function(){return this.getLabelPosition("y2",this.owner.config.axis_rotated?"inner-right":"inner-top")},f.prototype.getLabelPositionById=function(a){return"y2"===a?this.getY2AxisLabelPosition():"y"===a?this.getYAxisLabelPosition():this.getXAxisLabelPosition()},f.prototype.textForXAxisLabel=function(){return this.getLabelText("x")},f.prototype.textForYAxisLabel=function(){return this.getLabelText("y")},f.prototype.textForY2AxisLabel=function(){return this.getLabelText("y2")},f.prototype.xForAxisLabel=function(a,b){var c=this.owner;return a?b.isLeft?0:b.isCenter?c.width/2:c.width:b.isBottom?-c.height:b.isMiddle?-c.height/2:0},f.prototype.dxForAxisLabel=function(a,b){return a?b.isLeft?"0.5em":b.isRight?"-0.5em":"0":b.isTop?"-0.5em":b.isBottom?"0.5em":"0"},f.prototype.textAnchorForAxisLabel=function(a,b){return a?b.isLeft?"start":b.isCenter?"middle":"end":b.isBottom?"start":b.isMiddle?"middle":"end"},f.prototype.xForXAxisLabel=function(){return this.xForAxisLabel(!this.owner.config.axis_rotated,this.getXAxisLabelPosition())},f.prototype.xForYAxisLabel=function(){return this.xForAxisLabel(this.owner.config.axis_rotated,this.getYAxisLabelPosition())},f.prototype.xForY2AxisLabel=function(){return this.xForAxisLabel(this.owner.config.axis_rotated,this.getY2AxisLabelPosition())},f.prototype.dxForXAxisLabel=function(){return this.dxForAxisLabel(!this.owner.config.axis_rotated,this.getXAxisLabelPosition())},f.prototype.dxForYAxisLabel=function(){return this.dxForAxisLabel(this.owner.config.axis_rotated,this.getYAxisLabelPosition())},f.prototype.dxForY2AxisLabel=function(){return this.dxForAxisLabel(this.owner.config.axis_rotated,this.getY2AxisLabelPosition())},f.prototype.dyForXAxisLabel=function(){var a=this.owner,b=a.config,c=this.getXAxisLabelPosition();return b.axis_rotated?c.isInner?"1.2em":-25-this.getMaxTickWidth("x"):c.isInner?"-0.5em":b.axis_x_height?b.axis_x_height-10:"3em"},f.prototype.dyForYAxisLabel=function(){var a=this.owner,b=this.getYAxisLabelPosition();return a.config.axis_rotated?b.isInner?"-0.5em":"3em":b.isInner?"1.2em":-10-(a.config.axis_y_inner?0:this.getMaxTickWidth("y")+10)},f.prototype.dyForY2AxisLabel=function(){var a=this.owner,b=this.getY2AxisLabelPosition();return a.config.axis_rotated?b.isInner?"1.2em":"-2.2em":b.isInner?"-0.5em":15+(a.config.axis_y2_inner?0:this.getMaxTickWidth("y2")+15)},f.prototype.textAnchorForXAxisLabel=function(){var a=this.owner;return this.textAnchorForAxisLabel(!a.config.axis_rotated,this.getXAxisLabelPosition())},f.prototype.textAnchorForYAxisLabel=function(){var a=this.owner;return this.textAnchorForAxisLabel(a.config.axis_rotated,this.getYAxisLabelPosition())},f.prototype.textAnchorForY2AxisLabel=function(){var a=this.owner;return this.textAnchorForAxisLabel(a.config.axis_rotated,this.getY2AxisLabelPosition())},f.prototype.getMaxTickWidth=function(a,b){var c,d,e,f,g,h=this.owner,i=h.config,j=0;return b&&h.currentMaxTickWidths[a]?h.currentMaxTickWidths[a]:(h.svg&&(c=h.filterTargetsToShow(h.data.targets),"y"===a?(d=h.y.copy().domain(h.getYDomain(c,"y")),e=this.getYAxis(d,h.yOrient,i.axis_y_tick_format,h.yAxisTickValues,!1,!0)):"y2"===a?(d=h.y2.copy().domain(h.getYDomain(c,"y2")),e=this.getYAxis(d,h.y2Orient,i.axis_y2_tick_format,h.y2AxisTickValues,!1,!0)):(d=h.x.copy().domain(h.getXDomain(c)),e=this.getXAxis(d,h.xOrient,h.xAxisTickFormat,h.xAxisTickValues,!1,!0,!0),this.updateXAxisTickValues(c,e)),f=h.d3.select("body").append("div").classed("c3",!0),g=f.append("svg").style("visibility","hidden").style("position","fixed").style("top",0).style("left",0),g.append("g").call(e).each(function(){h.d3.select(this).selectAll("text").each(function(){var a=this.getBoundingClientRect();j<a.width&&(j=a.width)}),f.remove()})),h.currentMaxTickWidths[a]=0>=j?h.currentMaxTickWidths[a]:j,h.currentMaxTickWidths[a])},f.prototype.updateLabels=function(a){var b=this.owner,c=b.main.select("."+l.axisX+" ."+l.axisXLabel),d=b.main.select("."+l.axisY+" ."+l.axisYLabel),e=b.main.select("."+l.axisY2+" ."+l.axisY2Label);(a?c.transition():c).attr("x",this.xForXAxisLabel.bind(this)).attr("dx",this.dxForXAxisLabel.bind(this)).attr("dy",this.dyForXAxisLabel.bind(this)).text(this.textForXAxisLabel.bind(this)),(a?d.transition():d).attr("x",this.xForYAxisLabel.bind(this)).attr("dx",this.dxForYAxisLabel.bind(this)).attr("dy",this.dyForYAxisLabel.bind(this)).text(this.textForYAxisLabel.bind(this)),(a?e.transition():e).attr("x",this.xForY2AxisLabel.bind(this)).attr("dx",this.dxForY2AxisLabel.bind(this)).attr("dy",this.dyForY2AxisLabel.bind(this)).text(this.textForY2AxisLabel.bind(this))},f.prototype.getPadding=function(a,b,c,d){return m(a[b])?"ratio"===a.unit?a[b]*d:this.convertPixelsToAxisPadding(a[b],d):c},f.prototype.convertPixelsToAxisPadding=function(a,b){var c=this.owner,d=c.config.axis_rotated?c.width:c.height;return b*(a/d)},f.prototype.generateTickValues=function(a,b,c){var d,e,f,g,h,i,j,k=a;if(b)if(d=n(b)?b():b,1===d)k=[a[0]];else if(2===d)k=[a[0],a[a.length-1]];else if(d>2){for(g=d-2,e=a[0],f=a[a.length-1],h=(f-e)/(g+1),k=[e],i=0;g>i;i++)j=+e+h*(i+1),k.push(c?new Date(j):j);k.push(f)}return c||(k=k.sort(function(a,b){return a-b})),k},f.prototype.generateTransitions=function(a){var b=this.owner,c=b.axes;return{axisX:a?c.x.transition().duration(a):c.x,axisY:a?c.y.transition().duration(a):c.y,axisY2:a?c.y2.transition().duration(a):c.y2,axisSubX:a?c.subx.transition().duration(a):c.subx}},f.prototype.redraw=function(a,b){var c=this.owner;c.axes.x.style("opacity",b?0:1),c.axes.y.style("opacity",b?0:1),c.axes.y2.style("opacity",b?0:1),c.axes.subx.style("opacity",b?0:1),a.axisX.call(c.xAxis),a.axisY.call(c.yAxis),a.axisY2.call(c.y2Axis),a.axisSubX.call(c.subXAxis)},i.getClipPath=function(b){var c=a.navigator.appVersion.toLowerCase().indexOf("msie 9.")>=0;return"url("+(c?"":document.URL.split("#")[0])+"#"+b+")"},i.appendClip=function(a,b){return a.append("clipPath").attr("id",b).append("rect")},i.getAxisClipX=function(a){var b=Math.max(30,this.margin.left);return a?-(1+b):-(b-1)},i.getAxisClipY=function(a){return a?-20:-this.margin.top},i.getXAxisClipX=function(){var a=this;return a.getAxisClipX(!a.config.axis_rotated)},i.getXAxisClipY=function(){var a=this;return a.getAxisClipY(!a.config.axis_rotated)},i.getYAxisClipX=function(){var a=this;return a.config.axis_y_inner?-1:a.getAxisClipX(a.config.axis_rotated)},i.getYAxisClipY=function(){var a=this;return a.getAxisClipY(a.config.axis_rotated)},i.getAxisClipWidth=function(a){var b=this,c=Math.max(30,b.margin.left),d=Math.max(30,b.margin.right);return a?b.width+2+c+d:b.margin.left+20},i.getAxisClipHeight=function(a){return(a?this.margin.bottom:this.margin.top+this.height)+20},i.getXAxisClipWidth=function(){var a=this;return a.getAxisClipWidth(!a.config.axis_rotated)},i.getXAxisClipHeight=function(){var a=this;return a.getAxisClipHeight(!a.config.axis_rotated)},i.getYAxisClipWidth=function(){var a=this;return a.getAxisClipWidth(a.config.axis_rotated)+(a.config.axis_y_inner?20:0)},i.getYAxisClipHeight=function(){var a=this;return a.getAxisClipHeight(a.config.axis_rotated)},i.initPie=function(){var a=this,b=a.d3,c=a.config;a.pie=b.layout.pie().value(function(a){return a.values.reduce(function(a,b){return a+b.value},0)}),c.data_order||a.pie.sort(null)},i.updateRadius=function(){var a=this,b=a.config,c=b.gauge_width||b.donut_width;a.radiusExpanded=Math.min(a.arcWidth,a.arcHeight)/2,a.radius=.95*a.radiusExpanded,a.innerRadiusRatio=c?(a.radius-c)/a.radius:.6,a.innerRadius=a.hasType("donut")||a.hasType("gauge")?a.radius*a.innerRadiusRatio:0},i.updateArc=function(){var a=this;a.svgArc=a.getSvgArc(),a.svgArcExpanded=a.getSvgArcExpanded(),a.svgArcExpandedSub=a.getSvgArcExpanded(.98)},i.updateAngle=function(a){var b,c,d=this,e=d.config,f=!1,g=0,h=e.gauge_min,i=e.gauge_max;
return d.pie(d.filterTargetsToShow(d.data.targets)).forEach(function(b){f||b.data.id!==a.data.id||(f=!0,a=b,a.index=g),g++}),isNaN(a.startAngle)&&(a.startAngle=0),isNaN(a.endAngle)&&(a.endAngle=a.startAngle),d.isGaugeType(a.data)&&(b=Math.PI/(i-h),c=a.value<h?0:a.value<i?a.value-h:i-h,a.startAngle=-1*(Math.PI/2),a.endAngle=a.startAngle+b*c),f?a:null},i.getSvgArc=function(){var a=this,b=a.d3.svg.arc().outerRadius(a.radius).innerRadius(a.innerRadius),c=function(c,d){var e;return d?b(c):(e=a.updateAngle(c),e?b(e):"M 0 0")};return c.centroid=b.centroid,c},i.getSvgArcExpanded=function(a){var b=this,c=b.d3.svg.arc().outerRadius(b.radiusExpanded*(a?a:1)).innerRadius(b.innerRadius);return function(a){var d=b.updateAngle(a);return d?c(d):"M 0 0"}},i.getArc=function(a,b,c){return c||this.isArcType(a.data)?this.svgArc(a,b):"M 0 0"},i.transformForArcLabel=function(a){var b,c,d,e,f,g=this,h=g.updateAngle(a),i="";return h&&!g.hasType("gauge")&&(b=this.svgArc.centroid(h),c=isNaN(b[0])?0:b[0],d=isNaN(b[1])?0:b[1],e=Math.sqrt(c*c+d*d),f=g.radius&&e?(36/g.radius>.375?1.175-36/g.radius:.8)*g.radius/e:0,i="translate("+c*f+","+d*f+")"),i},i.getArcRatio=function(a){var b=this,c=b.hasType("gauge")?Math.PI:2*Math.PI;return a?(a.endAngle-a.startAngle)/c:null},i.convertToArcData=function(a){return this.addName({id:a.data.id,value:a.value,ratio:this.getArcRatio(a),index:a.index})},i.textForArcLabel=function(a){var b,c,d,e,f,g=this;return g.shouldShowArcLabel()?(b=g.updateAngle(a),c=b?b.value:null,d=g.getArcRatio(b),e=a.data.id,g.hasType("gauge")||g.meetsArcLabelThreshold(d)?(f=g.getArcLabelFormat(),f?f(c,d,e):g.defaultArcValueFormat(c,d)):""):""},i.expandArc=function(b){var c,d=this;return d.transiting?void(c=a.setInterval(function(){d.transiting||(a.clearInterval(c),d.legend.selectAll(".c3-legend-item-focused").size()>0&&d.expandArc(b))},10)):(b=d.mapToTargetIds(b),void d.svg.selectAll(d.selectorTargets(b,"."+l.chartArc)).each(function(a){d.shouldExpand(a.data.id)&&d.d3.select(this).selectAll("path").transition().duration(50).attr("d",d.svgArcExpanded).transition().duration(100).attr("d",d.svgArcExpandedSub).each(function(a){d.isDonutType(a.data)})}))},i.unexpandArc=function(a){var b=this;b.transiting||(a=b.mapToTargetIds(a),b.svg.selectAll(b.selectorTargets(a,"."+l.chartArc)).selectAll("path").transition().duration(50).attr("d",b.svgArc),b.svg.selectAll("."+l.arc).style("opacity",1))},i.shouldExpand=function(a){var b=this,c=b.config;return b.isDonutType(a)&&c.donut_expand||b.isGaugeType(a)&&c.gauge_expand||b.isPieType(a)&&c.pie_expand},i.shouldShowArcLabel=function(){var a=this,b=a.config,c=!0;return a.hasType("donut")?c=b.donut_label_show:a.hasType("pie")&&(c=b.pie_label_show),c},i.meetsArcLabelThreshold=function(a){var b=this,c=b.config,d=b.hasType("donut")?c.donut_label_threshold:c.pie_label_threshold;return a>=d},i.getArcLabelFormat=function(){var a=this,b=a.config,c=b.pie_label_format;return a.hasType("gauge")?c=b.gauge_label_format:a.hasType("donut")&&(c=b.donut_label_format),c},i.getArcTitle=function(){var a=this;return a.hasType("donut")?a.config.donut_title:""},i.updateTargetsForArc=function(a){var b,c,d=this,e=d.main,f=d.classChartArc.bind(d),g=d.classArcs.bind(d),h=d.classFocus.bind(d);b=e.select("."+l.chartArcs).selectAll("."+l.chartArc).data(d.pie(a)).attr("class",function(a){return f(a)+h(a.data)}),c=b.enter().append("g").attr("class",f),c.append("g").attr("class",g),c.append("text").attr("dy",d.hasType("gauge")?"-.1em":".35em").style("opacity",0).style("text-anchor","middle").style("pointer-events","none")},i.initArc=function(){var a=this;a.arcs=a.main.select("."+l.chart).append("g").attr("class",l.chartArcs).attr("transform",a.getTranslate("arc")),a.arcs.append("text").attr("class",l.chartArcsTitle).style("text-anchor","middle").text(a.getArcTitle())},i.redrawArc=function(a,b,c){var d,e=this,f=e.d3,g=e.config,h=e.main;d=h.selectAll("."+l.arcs).selectAll("."+l.arc).data(e.arcData.bind(e)),d.enter().append("path").attr("class",e.classArc.bind(e)).style("fill",function(a){return e.color(a.data)}).style("cursor",function(a){return g.interaction_enabled&&g.data_selection_isselectable(a)?"pointer":null}).style("opacity",0).each(function(a){e.isGaugeType(a.data)&&(a.startAngle=a.endAngle=-1*(Math.PI/2)),this._current=a}),d.attr("transform",function(a){return!e.isGaugeType(a.data)&&c?"scale(0)":""}).style("opacity",function(a){return a===this._current?0:1}).on("mouseover",g.interaction_enabled?function(a){var b,c;e.transiting||(b=e.updateAngle(a),c=e.convertToArcData(b),e.expandArc(b.data.id),e.api.focus(b.data.id),e.toggleFocusLegend(b.data.id,!0),e.config.data_onmouseover(c,this))}:null).on("mousemove",g.interaction_enabled?function(a){var b=e.updateAngle(a),c=e.convertToArcData(b),d=[c];e.showTooltip(d,this)}:null).on("mouseout",g.interaction_enabled?function(a){var b,c;e.transiting||(b=e.updateAngle(a),c=e.convertToArcData(b),e.unexpandArc(b.data.id),e.api.revert(),e.revertLegend(),e.hideTooltip(),e.config.data_onmouseout(c,this))}:null).on("click",g.interaction_enabled?function(a,b){var c=e.updateAngle(a),d=e.convertToArcData(c);e.toggleShape&&e.toggleShape(this,d,b),e.config.data_onclick.call(e.api,d,this)}:null).each(function(){e.transiting=!0}).transition().duration(a).attrTween("d",function(a){var b,c=e.updateAngle(a);return c?(isNaN(this._current.startAngle)&&(this._current.startAngle=0),isNaN(this._current.endAngle)&&(this._current.endAngle=this._current.startAngle),b=f.interpolate(this._current,c),this._current=b(0),function(c){var d=b(c);return d.data=a.data,e.getArc(d,!0)}):function(){return"M 0 0"}}).attr("transform",c?"scale(1)":"").style("fill",function(a){return e.levelColor?e.levelColor(a.data.values[0].value):e.color(a.data.id)}).style("opacity",1).call(e.endall,function(){e.transiting=!1}),d.exit().transition().duration(b).style("opacity",0).remove(),h.selectAll("."+l.chartArc).select("text").style("opacity",0).attr("class",function(a){return e.isGaugeType(a.data)?l.gaugeValue:""}).text(e.textForArcLabel.bind(e)).attr("transform",e.transformForArcLabel.bind(e)).style("font-size",function(a){return e.isGaugeType(a.data)?Math.round(e.radius/5)+"px":""}).transition().duration(a).style("opacity",function(a){return e.isTargetToShow(a.data.id)&&e.isArcType(a.data)?1:0}),h.select("."+l.chartArcsTitle).style("opacity",e.hasType("donut")||e.hasType("gauge")?1:0),e.hasType("gauge")&&(e.arcs.select("."+l.chartArcsBackground).attr("d",function(){var a={data:[{value:g.gauge_max}],startAngle:-1*(Math.PI/2),endAngle:Math.PI/2};return e.getArc(a,!0,!0)}),e.arcs.select("."+l.chartArcsGaugeUnit).attr("dy",".75em").text(g.gauge_label_show?g.gauge_units:""),e.arcs.select("."+l.chartArcsGaugeMin).attr("dx",-1*(e.innerRadius+(e.radius-e.innerRadius)/2)+"px").attr("dy","1.2em").text(g.gauge_label_show?g.gauge_min:""),e.arcs.select("."+l.chartArcsGaugeMax).attr("dx",e.innerRadius+(e.radius-e.innerRadius)/2+"px").attr("dy","1.2em").text(g.gauge_label_show?g.gauge_max:""))},i.initGauge=function(){var a=this.arcs;this.hasType("gauge")&&(a.append("path").attr("class",l.chartArcsBackground),a.append("text").attr("class",l.chartArcsGaugeUnit).style("text-anchor","middle").style("pointer-events","none"),a.append("text").attr("class",l.chartArcsGaugeMin).style("text-anchor","middle").style("pointer-events","none"),a.append("text").attr("class",l.chartArcsGaugeMax).style("text-anchor","middle").style("pointer-events","none"))},i.getGaugeLabelHeight=function(){return this.config.gauge_label_show?20:0},i.initRegion=function(){var a=this;a.region=a.main.append("g").attr("clip-path",a.clipPath).attr("class",l.regions)},i.updateRegion=function(a){var b=this,c=b.config;b.region.style("visibility",b.hasArcType()?"hidden":"visible"),b.mainRegion=b.main.select("."+l.regions).selectAll("."+l.region).data(c.regions),b.mainRegion.enter().append("g").attr("class",b.classRegion.bind(b)).append("rect").style("fill-opacity",0),b.mainRegion.exit().transition().duration(a).style("opacity",0).remove()},i.redrawRegion=function(a){var b=this,c=b.mainRegion.selectAll("rect"),d=b.regionX.bind(b),e=b.regionY.bind(b),f=b.regionWidth.bind(b),g=b.regionHeight.bind(b);return[(a?c.transition():c).attr("x",d).attr("y",e).attr("width",f).attr("height",g).style("fill-opacity",function(a){return m(a.opacity)?a.opacity:.1})]},i.regionX=function(a){var b,c=this,d=c.config,e="y"===a.axis?c.y:c.y2;return b="y"===a.axis||"y2"===a.axis?d.axis_rotated&&"start"in a?e(a.start):0:d.axis_rotated?0:"start"in a?c.x(c.isTimeSeries()?c.parseDate(a.start):a.start):0},i.regionY=function(a){var b,c=this,d=c.config,e="y"===a.axis?c.y:c.y2;return b="y"===a.axis||"y2"===a.axis?d.axis_rotated?0:"end"in a?e(a.end):0:d.axis_rotated&&"start"in a?c.x(c.isTimeSeries()?c.parseDate(a.start):a.start):0},i.regionWidth=function(a){var b,c=this,d=c.config,e=c.regionX(a),f="y"===a.axis?c.y:c.y2;return b="y"===a.axis||"y2"===a.axis?d.axis_rotated&&"end"in a?f(a.end):c.width:d.axis_rotated?c.width:"end"in a?c.x(c.isTimeSeries()?c.parseDate(a.end):a.end):c.width,e>b?0:b-e},i.regionHeight=function(a){var b,c=this,d=c.config,e=this.regionY(a),f="y"===a.axis?c.y:c.y2;return b="y"===a.axis||"y2"===a.axis?d.axis_rotated?c.height:"start"in a?f(a.start):c.height:d.axis_rotated&&"end"in a?c.x(c.isTimeSeries()?c.parseDate(a.end):a.end):c.height,e>b?0:b-e},i.isRegionOnX=function(a){return!a.axis||"x"===a.axis},i.drag=function(a){var b,c,d,e,f,g,h,i,j=this,k=j.config,m=j.main,n=j.d3;j.hasArcType()||k.data_selection_enabled&&(!k.zoom_enabled||j.zoom.altDomain)&&k.data_selection_multiple&&(b=j.dragStart[0],c=j.dragStart[1],d=a[0],e=a[1],f=Math.min(b,d),g=Math.max(b,d),h=k.data_selection_grouped?j.margin.top:Math.min(c,e),i=k.data_selection_grouped?j.height:Math.max(c,e),m.select("."+l.dragarea).attr("x",f).attr("y",h).attr("width",g-f).attr("height",i-h),m.selectAll("."+l.shapes).selectAll("."+l.shape).filter(function(a){return k.data_selection_isselectable(a)}).each(function(a,b){var c,d,e,k,m,o,p=n.select(this),q=p.classed(l.SELECTED),r=p.classed(l.INCLUDED),s=!1;if(p.classed(l.circle))c=1*p.attr("cx"),d=1*p.attr("cy"),m=j.togglePoint,s=c>f&&g>c&&d>h&&i>d;else{if(!p.classed(l.bar))return;o=y(this),c=o.x,d=o.y,e=o.width,k=o.height,m=j.togglePath,s=!(c>g||f>c+e||d>i||h>d+k)}s^r&&(p.classed(l.INCLUDED,!r),p.classed(l.SELECTED,!q),m.call(j,!q,p,a,b))}))},i.dragstart=function(a){var b=this,c=b.config;b.hasArcType()||c.data_selection_enabled&&(b.dragStart=a,b.main.select("."+l.chart).append("rect").attr("class",l.dragarea).style("opacity",.1),b.dragging=!0)},i.dragend=function(){var a=this,b=a.config;a.hasArcType()||b.data_selection_enabled&&(a.main.select("."+l.dragarea).transition().duration(100).style("opacity",0).remove(),a.main.selectAll("."+l.shape).classed(l.INCLUDED,!1),a.dragging=!1)},i.selectPoint=function(a,b,c){var d=this,e=d.config,f=(e.axis_rotated?d.circleY:d.circleX).bind(d),g=(e.axis_rotated?d.circleX:d.circleY).bind(d),h=d.pointSelectR.bind(d);e.data_onselected.call(d.api,b,a.node()),d.main.select("."+l.selectedCircles+d.getTargetSelectorSuffix(b.id)).selectAll("."+l.selectedCircle+"-"+c).data([b]).enter().append("circle").attr("class",function(){return d.generateClass(l.selectedCircle,c)}).attr("cx",f).attr("cy",g).attr("stroke",function(){return d.color(b)}).attr("r",function(a){return 1.4*d.pointSelectR(a)}).transition().duration(100).attr("r",h)},i.unselectPoint=function(a,b,c){var d=this;d.config.data_onunselected(b,a.node()),d.main.select("."+l.selectedCircles+d.getTargetSelectorSuffix(b.id)).selectAll("."+l.selectedCircle+"-"+c).transition().duration(100).attr("r",0).remove()},i.togglePoint=function(a,b,c,d){a?this.selectPoint(b,c,d):this.unselectPoint(b,c,d)},i.selectPath=function(a,b){var c=this;c.config.data_onselected.call(c,b,a.node()),a.transition().duration(100).style("fill",function(){return c.d3.rgb(c.color(b)).brighter(.75)})},i.unselectPath=function(a,b){var c=this;c.config.data_onunselected.call(c,b,a.node()),a.transition().duration(100).style("fill",function(){return c.color(b)})},i.togglePath=function(a,b,c,d){a?this.selectPath(b,c,d):this.unselectPath(b,c,d)},i.getToggle=function(a,b){var c,d=this;return"circle"===a.nodeName?c=d.isStepType(b)?function(){}:d.togglePoint:"path"===a.nodeName&&(c=d.togglePath),c},i.toggleShape=function(a,b,c){var d=this,e=d.d3,f=d.config,g=e.select(a),h=g.classed(l.SELECTED),i=d.getToggle(a,b).bind(d);f.data_selection_enabled&&f.data_selection_isselectable(b)&&(f.data_selection_multiple||d.main.selectAll("."+l.shapes+(f.data_selection_grouped?d.getTargetSelectorSuffix(b.id):"")).selectAll("."+l.shape).each(function(a,b){var c=e.select(this);c.classed(l.SELECTED)&&i(!1,c.classed(l.SELECTED,!1),a,b)}),g.classed(l.SELECTED,!h),i(!h,g,b,c))},i.initBrush=function(){var a=this,b=a.d3;a.brush=b.svg.brush().on("brush",function(){a.redrawForBrush()}),a.brush.update=function(){return a.context&&a.context.select("."+l.brush).call(this),this},a.brush.scale=function(b){return a.config.axis_rotated?this.y(b):this.x(b)}},i.initSubchart=function(){var a=this,b=a.config,c=a.context=a.svg.append("g").attr("transform",a.getTranslate("context"));c.style("visibility",b.subchart_show?"visible":"hidden"),c.append("g").attr("clip-path",a.clipPathForSubchart).attr("class",l.chart),c.select("."+l.chart).append("g").attr("class",l.chartBars),c.select("."+l.chart).append("g").attr("class",l.chartLines),c.append("g").attr("clip-path",a.clipPath).attr("class",l.brush).call(a.brush),a.axes.subx=c.append("g").attr("class",l.axisX).attr("transform",a.getTranslate("subx")).attr("clip-path",b.axis_rotated?"":a.clipPathForXAxis)},i.updateTargetsForSubchart=function(a){var b,c,d,e,f=this,g=f.context,h=f.config,i=f.classChartBar.bind(f),j=f.classBars.bind(f),k=f.classChartLine.bind(f),m=f.classLines.bind(f),n=f.classAreas.bind(f);h.subchart_show&&(e=g.select("."+l.chartBars).selectAll("."+l.chartBar).data(a).attr("class",i),d=e.enter().append("g").style("opacity",0).attr("class",i),d.append("g").attr("class",j),c=g.select("."+l.chartLines).selectAll("."+l.chartLine).data(a).attr("class",k),b=c.enter().append("g").style("opacity",0).attr("class",k),b.append("g").attr("class",m),b.append("g").attr("class",n),g.selectAll("."+l.brush+" rect").attr(h.axis_rotated?"width":"height",h.axis_rotated?f.width2:f.height2))},i.updateBarForSubchart=function(a){var b=this;b.contextBar=b.context.selectAll("."+l.bars).selectAll("."+l.bar).data(b.barData.bind(b)),b.contextBar.enter().append("path").attr("class",b.classBar.bind(b)).style("stroke","none").style("fill",b.color),b.contextBar.style("opacity",b.initialOpacity.bind(b)),b.contextBar.exit().transition().duration(a).style("opacity",0).remove()},i.redrawBarForSubchart=function(a,b,c){(b?this.contextBar.transition().duration(c):this.contextBar).attr("d",a).style("opacity",1)},i.updateLineForSubchart=function(a){var b=this;b.contextLine=b.context.selectAll("."+l.lines).selectAll("."+l.line).data(b.lineData.bind(b)),b.contextLine.enter().append("path").attr("class",b.classLine.bind(b)).style("stroke",b.color),b.contextLine.style("opacity",b.initialOpacity.bind(b)),b.contextLine.exit().transition().duration(a).style("opacity",0).remove()},i.redrawLineForSubchart=function(a,b,c){(b?this.contextLine.transition().duration(c):this.contextLine).attr("d",a).style("opacity",1)},i.updateAreaForSubchart=function(a){var b=this,c=b.d3;b.contextArea=b.context.selectAll("."+l.areas).selectAll("."+l.area).data(b.lineData.bind(b)),b.contextArea.enter().append("path").attr("class",b.classArea.bind(b)).style("fill",b.color).style("opacity",function(){return b.orgAreaOpacity=+c.select(this).style("opacity"),0}),b.contextArea.style("opacity",0),b.contextArea.exit().transition().duration(a).style("opacity",0).remove()},i.redrawAreaForSubchart=function(a,b,c){(b?this.contextArea.transition().duration(c):this.contextArea).attr("d",a).style("fill",this.color).style("opacity",this.orgAreaOpacity)},i.redrawSubchart=function(a,b,c,d,e,f,g){var h,i,j,k=this,l=k.d3,m=k.config;k.context.style("visibility",m.subchart_show?"visible":"hidden"),m.subchart_show&&(l.event&&"zoom"===l.event.type&&k.brush.extent(k.x.orgDomain()).update(),a&&(k.brush.empty()||k.brush.extent(k.x.orgDomain()).update(),h=k.generateDrawArea(e,!0),i=k.generateDrawBar(f,!0),j=k.generateDrawLine(g,!0),k.updateBarForSubchart(c),k.updateLineForSubchart(c),k.updateAreaForSubchart(c),k.redrawBarForSubchart(i,c,c),k.redrawLineForSubchart(j,c,c),k.redrawAreaForSubchart(h,c,c)))},i.redrawForBrush=function(){var a=this,b=a.x;a.redraw({withTransition:!1,withY:a.config.zoom_rescale,withSubchart:!1,withUpdateXDomain:!0,withDimension:!1}),a.config.subchart_onbrush.call(a.api,b.orgDomain())},i.transformContext=function(a,b){var c,d=this;b&&b.axisSubX?c=b.axisSubX:(c=d.context.select("."+l.axisX),a&&(c=c.transition())),d.context.attr("transform",d.getTranslate("context")),c.attr("transform",d.getTranslate("subx"))},i.getDefaultExtent=function(){var a=this,b=a.config,c=n(b.axis_x_extent)?b.axis_x_extent(a.getXDomain(a.data.targets)):b.axis_x_extent;return a.isTimeSeries()&&(c=[a.parseDate(c[0]),a.parseDate(c[1])]),c},i.initZoom=function(){var a,b=this,c=b.d3,d=b.config;b.zoom=c.behavior.zoom().on("zoomstart",function(){a=c.event.sourceEvent,b.zoom.altDomain=c.event.sourceEvent.altKey?b.x.orgDomain():null,d.zoom_onzoomstart.call(b.api,c.event.sourceEvent)}).on("zoom",function(){b.redrawForZoom.call(b)}).on("zoomend",function(){var e=c.event.sourceEvent;e&&a.clientX===e.clientX&&a.clientY===e.clientY||(b.redrawEventRect(),b.updateZoom(),d.zoom_onzoomend.call(b.api,b.x.orgDomain()))}),b.zoom.scale=function(a){return d.axis_rotated?this.y(a):this.x(a)},b.zoom.orgScaleExtent=function(){var a=d.zoom_extent?d.zoom_extent:[1,10];return[a[0],Math.max(b.getMaxDataCount()/a[1],a[1])]},b.zoom.updateScaleExtent=function(){var a=t(b.x.orgDomain())/t(b.orgXDomain),c=this.orgScaleExtent();return this.scaleExtent([c[0]*a,c[1]*a]),this}},i.updateZoom=function(){var a=this,b=a.config.zoom_enabled?a.zoom:function(){};a.main.select("."+l.zoomRect).call(b).on("dblclick.zoom",null),a.main.selectAll("."+l.eventRect).call(b).on("dblclick.zoom",null)},i.redrawForZoom=function(){var a=this,b=a.d3,c=a.config,d=a.zoom,e=a.x;if(c.zoom_enabled&&0!==a.filterTargetsToShow(a.data.targets).length){if("mousemove"===b.event.sourceEvent.type&&d.altDomain)return e.domain(d.altDomain),void d.scale(e).updateScaleExtent();a.isCategorized()&&e.orgDomain()[0]===a.orgXDomain[0]&&e.domain([a.orgXDomain[0]-1e-10,e.orgDomain()[1]]),a.redraw({withTransition:!1,withY:c.zoom_rescale,withSubchart:!1,withEventRect:!1,withDimension:!1}),"mousemove"===b.event.sourceEvent.type&&(a.cancelClick=!0),c.zoom_onzoom.call(a.api,e.orgDomain())}},i.generateColor=function(){var a=this,b=a.config,c=a.d3,d=b.data_colors,e=v(b.color_pattern)?b.color_pattern:c.scale.category10().range(),f=b.data_color,g=[];return function(a){var b,c=a.id||a.data&&a.data.id||a;return d[c]instanceof Function?b=d[c](a):d[c]?b=d[c]:(g.indexOf(c)<0&&g.push(c),b=e[g.indexOf(c)%e.length],d[c]=b),f instanceof Function?f(b,a):b}},i.generateLevelColor=function(){var a=this,b=a.config,c=b.color_pattern,d=b.color_threshold,e="value"===d.unit,f=d.values&&d.values.length?d.values:[],g=d.max||100;return v(b.color_threshold)?function(a){var b,d,h=c[c.length-1];for(b=0;b<f.length;b++)if(d=e?a:100*a/g,d<f[b]){h=c[b];break}return h}:null},i.getYFormat=function(a){var b=this,c=a&&!b.hasType("gauge")?b.defaultArcValueFormat:b.yFormat,d=a&&!b.hasType("gauge")?b.defaultArcValueFormat:b.y2Format;return function(a,e,f){var g="y2"===b.axis.getId(f)?d:c;return g.call(b,a,e)}},i.yFormat=function(a){var b=this,c=b.config,d=c.axis_y_tick_format?c.axis_y_tick_format:b.defaultValueFormat;return d(a)},i.y2Format=function(a){var b=this,c=b.config,d=c.axis_y2_tick_format?c.axis_y2_tick_format:b.defaultValueFormat;return d(a)},i.defaultValueFormat=function(a){return m(a)?+a:""},i.defaultArcValueFormat=function(a,b){return(100*b).toFixed(1)+"%"},i.dataLabelFormat=function(a){var b,c=this,d=c.config.data_labels,e=function(a){return m(a)?+a:""};return b="function"==typeof d.format?d.format:"object"==typeof d.format?d.format[a]?d.format[a]===!0?e:d.format[a]:function(){return""}:e},i.hasCaches=function(a){for(var b=0;b<a.length;b++)if(!(a[b]in this.cache))return!1;return!0},i.addCache=function(a,b){this.cache[a]=this.cloneTarget(b)},i.getCaches=function(a){var b,c=[];for(b=0;b<a.length;b++)a[b]in this.cache&&c.push(this.cloneTarget(this.cache[a[b]]));return c};var l=i.CLASS={target:"c3-target",chart:"c3-chart",chartLine:"c3-chart-line",chartLines:"c3-chart-lines",chartBar:"c3-chart-bar",chartBars:"c3-chart-bars",chartText:"c3-chart-text",chartTexts:"c3-chart-texts",chartArc:"c3-chart-arc",chartArcs:"c3-chart-arcs",chartArcsTitle:"c3-chart-arcs-title",chartArcsBackground:"c3-chart-arcs-background",chartArcsGaugeUnit:"c3-chart-arcs-gauge-unit",chartArcsGaugeMax:"c3-chart-arcs-gauge-max",chartArcsGaugeMin:"c3-chart-arcs-gauge-min",selectedCircle:"c3-selected-circle",selectedCircles:"c3-selected-circles",eventRect:"c3-event-rect",eventRects:"c3-event-rects",eventRectsSingle:"c3-event-rects-single",eventRectsMultiple:"c3-event-rects-multiple",zoomRect:"c3-zoom-rect",brush:"c3-brush",focused:"c3-focused",defocused:"c3-defocused",region:"c3-region",regions:"c3-regions",tooltipContainer:"c3-tooltip-container",tooltip:"c3-tooltip",tooltipName:"c3-tooltip-name",shape:"c3-shape",shapes:"c3-shapes",line:"c3-line",lines:"c3-lines",bar:"c3-bar",bars:"c3-bars",circle:"c3-circle",circles:"c3-circles",arc:"c3-arc",arcs:"c3-arcs",area:"c3-area",areas:"c3-areas",empty:"c3-empty",text:"c3-text",texts:"c3-texts",gaugeValue:"c3-gauge-value",grid:"c3-grid",gridLines:"c3-grid-lines",xgrid:"c3-xgrid",xgrids:"c3-xgrids",xgridLine:"c3-xgrid-line",xgridLines:"c3-xgrid-lines",xgridFocus:"c3-xgrid-focus",ygrid:"c3-ygrid",ygrids:"c3-ygrids",ygridLine:"c3-ygrid-line",ygridLines:"c3-ygrid-lines",axis:"c3-axis",axisX:"c3-axis-x",axisXLabel:"c3-axis-x-label",axisY:"c3-axis-y",axisYLabel:"c3-axis-y-label",axisY2:"c3-axis-y2",axisY2Label:"c3-axis-y2-label",legendBackground:"c3-legend-background",legendItem:"c3-legend-item",legendItemEvent:"c3-legend-item-event",legendItemTile:"c3-legend-item-tile",legendItemHidden:"c3-legend-item-hidden",legendItemFocused:"c3-legend-item-focused",dragarea:"c3-dragarea",EXPANDED:"_expanded_",SELECTED:"_selected_",INCLUDED:"_included_"};i.generateClass=function(a,b){return" "+a+" "+a+this.getTargetSelectorSuffix(b)},i.classText=function(a){return this.generateClass(l.text,a.index)},i.classTexts=function(a){return this.generateClass(l.texts,a.id)},i.classShape=function(a){return this.generateClass(l.shape,a.index)},i.classShapes=function(a){return this.generateClass(l.shapes,a.id)},i.classLine=function(a){return this.classShape(a)+this.generateClass(l.line,a.id)},i.classLines=function(a){return this.classShapes(a)+this.generateClass(l.lines,a.id)},i.classCircle=function(a){return this.classShape(a)+this.generateClass(l.circle,a.index)},i.classCircles=function(a){return this.classShapes(a)+this.generateClass(l.circles,a.id)},i.classBar=function(a){return this.classShape(a)+this.generateClass(l.bar,a.index)},i.classBars=function(a){return this.classShapes(a)+this.generateClass(l.bars,a.id)},i.classArc=function(a){return this.classShape(a.data)+this.generateClass(l.arc,a.data.id)},i.classArcs=function(a){return this.classShapes(a.data)+this.generateClass(l.arcs,a.data.id)},i.classArea=function(a){return this.classShape(a)+this.generateClass(l.area,a.id)},i.classAreas=function(a){return this.classShapes(a)+this.generateClass(l.areas,a.id)},i.classRegion=function(a,b){return this.generateClass(l.region,b)+" "+("class"in a?a["class"]:"")},i.classEvent=function(a){return this.generateClass(l.eventRect,a.index)},i.classTarget=function(a){var b=this,c=b.config.data_classes[a],d="";return c&&(d=" "+l.target+"-"+c),b.generateClass(l.target,a)+d},i.classFocus=function(a){return this.classFocused(a)+this.classDefocused(a)},i.classFocused=function(a){return" "+(this.focusedTargetIds.indexOf(a.id)>=0?l.focused:"")},i.classDefocused=function(a){return" "+(this.defocusedTargetIds.indexOf(a.id)>=0?l.defocused:"")},i.classChartText=function(a){return l.chartText+this.classTarget(a.id)},i.classChartLine=function(a){return l.chartLine+this.classTarget(a.id)},i.classChartBar=function(a){return l.chartBar+this.classTarget(a.id)},i.classChartArc=function(a){return l.chartArc+this.classTarget(a.data.id)},i.getTargetSelectorSuffix=function(a){return a||0===a?("-"+a).replace(/[\s?!@#$%^&*()_=+,.<>'":;\[\]\/|~`{}\\]/g,"-"):""},i.selectorTarget=function(a,b){return(b||"")+"."+l.target+this.getTargetSelectorSuffix(a)},i.selectorTargets=function(a,b){var c=this;return a=a||[],a.length?a.map(function(a){return c.selectorTarget(a,b)}):null},i.selectorLegend=function(a){return"."+l.legendItem+this.getTargetSelectorSuffix(a)},i.selectorLegends=function(a){var b=this;return a&&a.length?a.map(function(a){return b.selectorLegend(a)}):null};var m=i.isValue=function(a){return a||0===a},n=i.isFunction=function(a){return"function"==typeof a},o=i.isString=function(a){return"string"==typeof a},p=i.isUndefined=function(a){return"undefined"==typeof a},q=i.isDefined=function(a){return"undefined"!=typeof a},r=i.ceil10=function(a){return 10*Math.ceil(a/10)},s=i.asHalfPixel=function(a){return Math.ceil(a)+.5},t=i.diffDomain=function(a){return a[1]-a[0]},u=i.isEmpty=function(a){return!a||o(a)&&0===a.length||"object"==typeof a&&0===Object.keys(a).length},v=i.notEmpty=function(a){return Object.keys(a).length>0},w=i.getOption=function(a,b,c){return q(a[b])?a[b]:c},x=i.hasValue=function(a,b){var c=!1;return Object.keys(a).forEach(function(d){a[d]===b&&(c=!0)}),c},y=i.getPathBox=function(a){var b=a.getBoundingClientRect(),c=[a.pathSegList.getItem(0),a.pathSegList.getItem(1)],d=c[0].x,e=Math.min(c[0].y,c[1].y);return{x:d,y:e,width:b.width,height:b.height}};h.focus=function(a){var b,c=this.internal;a=c.mapToTargetIds(a),b=c.svg.selectAll(c.selectorTargets(a.filter(c.isTargetToShow,c))),this.revert(),this.defocus(),b.classed(l.focused,!0).classed(l.defocused,!1),c.hasArcType()&&c.expandArc(a),c.toggleFocusLegend(a,!0),c.focusedTargetIds=a,c.defocusedTargetIds=c.defocusedTargetIds.filter(function(b){return a.indexOf(b)<0})},h.defocus=function(a){var b,c=this.internal;a=c.mapToTargetIds(a),b=c.svg.selectAll(c.selectorTargets(a.filter(c.isTargetToShow,c))),b.classed(l.focused,!1).classed(l.defocused,!0),c.hasArcType()&&c.unexpandArc(a),c.toggleFocusLegend(a,!1),c.focusedTargetIds=c.focusedTargetIds.filter(function(b){return a.indexOf(b)<0}),c.defocusedTargetIds=a},h.revert=function(a){var b,c=this.internal;a=c.mapToTargetIds(a),b=c.svg.selectAll(c.selectorTargets(a)),b.classed(l.focused,!1).classed(l.defocused,!1),c.hasArcType()&&c.unexpandArc(a),c.config.legend_show&&(c.showLegend(a.filter(c.isLegendToShow.bind(c))),c.legend.selectAll(c.selectorLegends(a)).filter(function(){return c.d3.select(this).classed(l.legendItemFocused)}).classed(l.legendItemFocused,!1)),c.focusedTargetIds=[],c.defocusedTargetIds=[]},h.show=function(a,b){var c,d=this.internal;a=d.mapToTargetIds(a),b=b||{},d.removeHiddenTargetIds(a),c=d.svg.selectAll(d.selectorTargets(a)),c.transition().style("opacity",1,"important").call(d.endall,function(){c.style("opacity",null).style("opacity",1)}),b.withLegend&&d.showLegend(a),d.redraw({withUpdateOrgXDomain:!0,withUpdateXDomain:!0,withLegend:!0})},h.hide=function(a,b){var c,d=this.internal;a=d.mapToTargetIds(a),b=b||{},d.addHiddenTargetIds(a),c=d.svg.selectAll(d.selectorTargets(a)),c.transition().style("opacity",0,"important").call(d.endall,function(){c.style("opacity",null).style("opacity",0)}),b.withLegend&&d.hideLegend(a),d.redraw({withUpdateOrgXDomain:!0,withUpdateXDomain:!0,withLegend:!0})},h.toggle=function(a,b){var c=this,d=this.internal;d.mapToTargetIds(a).forEach(function(a){d.isTargetToShow(a)?c.hide(a,b):c.show(a,b)})},h.zoom=function(a){var b=this.internal;return a&&(b.isTimeSeries()&&(a=a.map(function(a){return b.parseDate(a)})),b.brush.extent(a),b.redraw({withUpdateXDomain:!0,withY:b.config.zoom_rescale}),b.config.zoom_onzoom.call(this,b.x.orgDomain())),b.brush.extent()},h.zoom.enable=function(a){var b=this.internal;b.config.zoom_enabled=a,b.updateAndRedraw()},h.unzoom=function(){var a=this.internal;a.brush.clear().update(),a.redraw({withUpdateXDomain:!0})},h.load=function(a){var b=this.internal,c=b.config;return a.xs&&b.addXs(a.xs),"classes"in a&&Object.keys(a.classes).forEach(function(b){c.data_classes[b]=a.classes[b]}),"categories"in a&&b.isCategorized()&&(c.axis_x_categories=a.categories),"axes"in a&&Object.keys(a.axes).forEach(function(b){c.data_axes[b]=a.axes[b]}),"colors"in a&&Object.keys(a.colors).forEach(function(b){c.data_colors[b]=a.colors[b]}),"cacheIds"in a&&b.hasCaches(a.cacheIds)?void b.load(b.getCaches(a.cacheIds),a.done):void("unload"in a?b.unload(b.mapToTargetIds("boolean"==typeof a.unload&&a.unload?null:a.unload),function(){b.loadFromArgs(a)}):b.loadFromArgs(a))},h.unload=function(a){var b=this.internal;a=a||{},a instanceof Array?a={ids:a}:"string"==typeof a&&(a={ids:[a]}),b.unload(b.mapToTargetIds(a.ids),function(){b.redraw({withUpdateOrgXDomain:!0,withUpdateXDomain:!0,withLegend:!0}),a.done&&a.done()})},h.flow=function(a){var b,c,d,e,f,g,h,i,j=this.internal,k=[],l=j.getMaxDataCount(),n=0,o=0;if(a.json)c=j.convertJsonToData(a.json,a.keys);else if(a.rows)c=j.convertRowsToData(a.rows);else{if(!a.columns)return;c=j.convertColumnsToData(a.columns)}b=j.convertDataToTargets(c,!0),j.data.targets.forEach(function(a){var c,d,e=!1;for(c=0;c<b.length;c++)if(a.id===b[c].id){for(e=!0,a.values[a.values.length-1]&&(o=a.values[a.values.length-1].index+1),n=b[c].values.length,d=0;n>d;d++)b[c].values[d].index=o+d,j.isTimeSeries()||(b[c].values[d].x=o+d);a.values=a.values.concat(b[c].values),b.splice(c,1);break}e||k.push(a.id)}),j.data.targets.forEach(function(a){var b,c;for(b=0;b<k.length;b++)if(a.id===k[b])for(o=a.values[a.values.length-1].index+1,c=0;n>c;c++)a.values.push({id:a.id,index:o+c,x:j.isTimeSeries()?j.getOtherTargetX(o+c):o+c,value:null})}),j.data.targets.length&&b.forEach(function(a){var b,c=[];for(b=j.data.targets[0].values[0].index;o>b;b++)c.push({id:a.id,index:b,x:j.isTimeSeries()?j.getOtherTargetX(b):b,value:null});a.values.forEach(function(a){a.index+=o,j.isTimeSeries()||(a.x+=o)}),a.values=c.concat(a.values)}),j.data.targets=j.data.targets.concat(b),d=j.getMaxDataCount(),f=j.data.targets[0],g=f.values[0],q(a.to)?(n=0,i=j.isTimeSeries()?j.parseDate(a.to):a.to,f.values.forEach(function(a){a.x<i&&n++})):q(a.length)&&(n=a.length),l?1===l&&j.isTimeSeries()&&(h=(f.values[f.values.length-1].x-g.x)/2,e=[new Date(+g.x-h),new Date(+g.x+h)],j.updateXDomain(null,!0,!0,!1,e)):(h=j.isTimeSeries()?f.values.length>1?f.values[f.values.length-1].x-g.x:g.x-j.getXDomain(j.data.targets)[0]:1,e=[g.x-h,g.x],j.updateXDomain(null,!0,!0,!1,e)),j.updateTargets(j.data.targets),j.redraw({flow:{index:g.index,length:n,duration:m(a.duration)?a.duration:j.config.transition_duration,done:a.done,orgDataCount:l},withLegend:!0,withTransition:l>1,withTrimXDomain:!1,withUpdateXAxis:!0})},i.generateFlow=function(a){var b=this,c=b.config,d=b.d3;return function(){var e,f,g,h=a.targets,i=a.flow,j=a.drawBar,k=a.drawLine,m=a.drawArea,n=a.cx,o=a.cy,p=a.xv,q=a.xForText,r=a.yForText,s=a.duration,u=1,v=i.index,w=i.length,x=b.getValueOnIndex(b.data.targets[0].values,v),y=b.getValueOnIndex(b.data.targets[0].values,v+w),z=b.x.domain(),A=i.duration||s,B=i.done||function(){},C=b.generateWait(),D=b.xgrid||d.selectAll([]),E=b.xgridLines||d.selectAll([]),F=b.mainRegion||d.selectAll([]),G=b.mainText||d.selectAll([]),H=b.mainBar||d.selectAll([]),I=b.mainLine||d.selectAll([]),J=b.mainArea||d.selectAll([]),K=b.mainCircle||d.selectAll([]);b.flowing=!0,b.data.targets.forEach(function(a){a.values.splice(0,w)}),g=b.updateXDomain(h,!0,!0),b.updateXGrid&&b.updateXGrid(!0),i.orgDataCount?e=1===i.orgDataCount||x.x===y.x?b.x(z[0])-b.x(g[0]):b.isTimeSeries()?b.x(z[0])-b.x(g[0]):b.x(x.x)-b.x(y.x):1!==b.data.targets[0].values.length?e=b.x(z[0])-b.x(g[0]):b.isTimeSeries()?(x=b.getValueOnIndex(b.data.targets[0].values,0),y=b.getValueOnIndex(b.data.targets[0].values,b.data.targets[0].values.length-1),e=b.x(x.x)-b.x(y.x)):e=t(g)/2,u=t(z)/t(g),f="translate("+e+",0) scale("+u+",1)",b.hideXGridFocus(),b.hideTooltip(),d.transition().ease("linear").duration(A).each(function(){C.add(b.axes.x.transition().call(b.xAxis)),C.add(H.transition().attr("transform",f)),C.add(I.transition().attr("transform",f)),C.add(J.transition().attr("transform",f)),C.add(K.transition().attr("transform",f)),C.add(G.transition().attr("transform",f)),C.add(F.filter(b.isRegionOnX).transition().attr("transform",f)),C.add(D.transition().attr("transform",f)),C.add(E.transition().attr("transform",f))
}).call(C,function(){var a,d=[],e=[],f=[];if(w){for(a=0;w>a;a++)d.push("."+l.shape+"-"+(v+a)),e.push("."+l.text+"-"+(v+a)),f.push("."+l.eventRect+"-"+(v+a));b.svg.selectAll("."+l.shapes).selectAll(d).remove(),b.svg.selectAll("."+l.texts).selectAll(e).remove(),b.svg.selectAll("."+l.eventRects).selectAll(f).remove(),b.svg.select("."+l.xgrid).remove()}D.attr("transform",null).attr(b.xgridAttr),E.attr("transform",null),E.select("line").attr("x1",c.axis_rotated?0:p).attr("x2",c.axis_rotated?b.width:p),E.select("text").attr("x",c.axis_rotated?b.width:0).attr("y",p),H.attr("transform",null).attr("d",j),I.attr("transform",null).attr("d",k),J.attr("transform",null).attr("d",m),K.attr("transform",null).attr("cx",n).attr("cy",o),G.attr("transform",null).attr("x",q).attr("y",r).style("fill-opacity",b.opacityForText.bind(b)),F.attr("transform",null),F.select("rect").filter(b.isRegionOnX).attr("x",b.regionX.bind(b)).attr("width",b.regionWidth.bind(b)),c.interaction_enabled&&b.redrawEventRect(),B(),b.flowing=!1})}},h.selected=function(a){var b=this.internal,c=b.d3;return c.merge(b.main.selectAll("."+l.shapes+b.getTargetSelectorSuffix(a)).selectAll("."+l.shape).filter(function(){return c.select(this).classed(l.SELECTED)}).map(function(a){return a.map(function(a){var b=a.__data__;return b.data?b.data:b})}))},h.select=function(a,b,c){var d=this.internal,e=d.d3,f=d.config;f.data_selection_enabled&&d.main.selectAll("."+l.shapes).selectAll("."+l.shape).each(function(g,h){var i=e.select(this),j=g.data?g.data.id:g.id,k=d.getToggle(this,g).bind(d),m=f.data_selection_grouped||!a||a.indexOf(j)>=0,n=!b||b.indexOf(h)>=0,o=i.classed(l.SELECTED);i.classed(l.line)||i.classed(l.area)||(m&&n?f.data_selection_isselectable(g)&&!o&&k(!0,i.classed(l.SELECTED,!0),g,h):q(c)&&c&&o&&k(!1,i.classed(l.SELECTED,!1),g,h))})},h.unselect=function(a,b){var c=this.internal,d=c.d3,e=c.config;e.data_selection_enabled&&c.main.selectAll("."+l.shapes).selectAll("."+l.shape).each(function(f,g){var h=d.select(this),i=f.data?f.data.id:f.id,j=c.getToggle(this,f).bind(c),k=e.data_selection_grouped||!a||a.indexOf(i)>=0,m=!b||b.indexOf(g)>=0,n=h.classed(l.SELECTED);h.classed(l.line)||h.classed(l.area)||k&&m&&e.data_selection_isselectable(f)&&n&&j(!1,h.classed(l.SELECTED,!1),f,g)})},h.transform=function(a,b){var c=this.internal,d=["pie","donut"].indexOf(a)>=0?{withTransform:!0}:null;c.transformTo(b,a,d)},i.transformTo=function(a,b,c){var d=this,e=!d.hasArcType(),f=c||{withTransitionForAxis:e};f.withTransitionForTransform=!1,d.transiting=!1,d.setTargetType(a,b),d.updateTargets(d.data.targets),d.updateAndRedraw(f)},h.groups=function(a){var b=this.internal,c=b.config;return p(a)?c.data_groups:(c.data_groups=a,b.redraw(),c.data_groups)},h.xgrids=function(a){var b=this.internal,c=b.config;return a?(c.grid_x_lines=a,b.redrawWithoutRescale(),c.grid_x_lines):c.grid_x_lines},h.xgrids.add=function(a){var b=this.internal;return this.xgrids(b.config.grid_x_lines.concat(a?a:[]))},h.xgrids.remove=function(a){var b=this.internal;b.removeGridLines(a,!0)},h.ygrids=function(a){var b=this.internal,c=b.config;return a?(c.grid_y_lines=a,b.redrawWithoutRescale(),c.grid_y_lines):c.grid_y_lines},h.ygrids.add=function(a){var b=this.internal;return this.ygrids(b.config.grid_y_lines.concat(a?a:[]))},h.ygrids.remove=function(a){var b=this.internal;b.removeGridLines(a,!1)},h.regions=function(a){var b=this.internal,c=b.config;return a?(c.regions=a,b.redrawWithoutRescale(),c.regions):c.regions},h.regions.add=function(a){var b=this.internal,c=b.config;return a?(c.regions=c.regions.concat(a),b.redrawWithoutRescale(),c.regions):c.regions},h.regions.remove=function(a){var b,c,d,e=this.internal,f=e.config;return a=a||{},b=e.getOption(a,"duration",f.transition_duration),c=e.getOption(a,"classes",[l.region]),d=e.main.select("."+l.regions).selectAll(c.map(function(a){return"."+a})),(b?d.transition().duration(b):d).style("opacity",0).remove(),f.regions=f.regions.filter(function(a){var b=!1;return a["class"]?(a["class"].split(" ").forEach(function(a){c.indexOf(a)>=0&&(b=!0)}),!b):!0}),f.regions},h.data=function(a){var b=this.internal.data.targets;return"undefined"==typeof a?b:b.filter(function(b){return[].concat(a).indexOf(b.id)>=0})},h.data.shown=function(a){return this.internal.filterTargetsToShow(this.data(a))},h.data.values=function(a){var b,c=null;return a&&(b=this.data(a),c=b[0]?b[0].values.map(function(a){return a.value}):null),c},h.data.names=function(a){return this.internal.clearLegendItemTextBoxCache(),this.internal.updateDataAttributes("names",a)},h.data.colors=function(a){return this.internal.updateDataAttributes("colors",a)},h.data.axes=function(a){return this.internal.updateDataAttributes("axes",a)},h.category=function(a,b){var c=this.internal,d=c.config;return arguments.length>1&&(d.axis_x_categories[a]=b,c.redraw()),d.axis_x_categories[a]},h.categories=function(a){var b=this.internal,c=b.config;return arguments.length?(c.axis_x_categories=a,b.redraw(),c.axis_x_categories):c.axis_x_categories},h.color=function(a){var b=this.internal;return b.color(a)},h.x=function(a){var b=this.internal;return arguments.length&&(b.updateTargetX(b.data.targets,a),b.redraw({withUpdateOrgXDomain:!0,withUpdateXDomain:!0})),b.data.xs},h.xs=function(a){var b=this.internal;return arguments.length&&(b.updateTargetXs(b.data.targets,a),b.redraw({withUpdateOrgXDomain:!0,withUpdateXDomain:!0})),b.data.xs},h.axis=function(){},h.axis.labels=function(a){var b=this.internal;arguments.length&&(Object.keys(a).forEach(function(c){b.axis.setLabelText(c,a[c])}),b.axis.updateLabels())},h.axis.max=function(a){var b=this.internal,c=b.config;return arguments.length?("object"==typeof a?(m(a.x)&&(c.axis_x_max=a.x),m(a.y)&&(c.axis_y_max=a.y),m(a.y2)&&(c.axis_y2_max=a.y2)):c.axis_y_max=c.axis_y2_max=a,void b.redraw({withUpdateOrgXDomain:!0,withUpdateXDomain:!0})):{x:c.axis_x_max,y:c.axis_y_max,y2:c.axis_y2_max}},h.axis.min=function(a){var b=this.internal,c=b.config;return arguments.length?("object"==typeof a?(m(a.x)&&(c.axis_x_min=a.x),m(a.y)&&(c.axis_y_min=a.y),m(a.y2)&&(c.axis_y2_min=a.y2)):c.axis_y_min=c.axis_y2_min=a,void b.redraw({withUpdateOrgXDomain:!0,withUpdateXDomain:!0})):{x:c.axis_x_min,y:c.axis_y_min,y2:c.axis_y2_min}},h.axis.range=function(a){return arguments.length?(q(a.max)&&this.axis.max(a.max),void(q(a.min)&&this.axis.min(a.min))):{max:this.axis.max(),min:this.axis.min()}},h.legend=function(){},h.legend.show=function(a){var b=this.internal;b.showLegend(b.mapToTargetIds(a)),b.updateAndRedraw({withLegend:!0})},h.legend.hide=function(a){var b=this.internal;b.hideLegend(b.mapToTargetIds(a)),b.updateAndRedraw({withLegend:!0})},h.resize=function(a){var b=this.internal,c=b.config;c.size_width=a?a.width:null,c.size_height=a?a.height:null,this.flush()},h.flush=function(){var a=this.internal;a.updateAndRedraw({withLegend:!0,withTransition:!1,withTransitionForTransform:!1})},h.destroy=function(){var b=this.internal;return a.clearInterval(b.intervalForObserveInserted),a.onresize=null,b.selectChart.classed("c3",!1).html(""),Object.keys(b).forEach(function(a){b[a]=null}),null},h.tooltip=function(){},h.tooltip.show=function(a){var b,c,d=this.internal;a.mouse&&(c=a.mouse),a.data?d.isMultipleX()?(c=[d.x(a.data.x),d.getYScale(a.data.id)(a.data.value)],b=null):b=m(a.data.index)?a.data.index:d.getIndexByX(a.data.x):"undefined"!=typeof a.x?b=d.getIndexByX(a.x):"undefined"!=typeof a.index&&(b=a.index),d.dispatchEvent("mouseover",b,c),d.dispatchEvent("mousemove",b,c)},h.tooltip.hide=function(){this.internal.dispatchEvent("mouseout",0)};var z;i.isSafari=function(){var b=a.navigator.userAgent;return b.indexOf("Safari")>=0&&b.indexOf("Chrome")<0},i.isChrome=function(){var b=a.navigator.userAgent;return b.indexOf("Chrome")>=0},Function.prototype.bind||(Function.prototype.bind=function(a){if("function"!=typeof this)throw new TypeError("Function.prototype.bind - what is trying to be bound is not callable");var b=Array.prototype.slice.call(arguments,1),c=this,d=function(){},e=function(){return c.apply(this instanceof d?this:a,b.concat(Array.prototype.slice.call(arguments)))};return d.prototype=this.prototype,e.prototype=new d,e}),"function"==typeof define&&define.amd?define("c3",["d3"],k):"undefined"!=typeof exports&&"undefined"!=typeof module?module.exports=k:a.c3=k}(window);
/**
 Copyright (c) 2014 BrightPoint Consulting, Inc.

 Permission is hereby granted, free of charge, to any person
 obtaining a copy of this software and associated documentation
 files (the "Software"), to deal in the Software without
 restriction, including without limitation the rights to use,
 copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following
 conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 OTHER DEALINGS IN THE SOFTWARE.
 */

function autocomplete(parent) {
    var _data=null,
        _delay= 0,
        _selection,
        _margin = {top: 30, right: 10, bottom: 50, left: 80},
        __width = 420,
        __height = 420,
        _placeHolder = "Search",
        _width,
        _height,
        _matches,
        _searchTerm,
        _lastSearchTerm,
        _currentIndex,
        _keys,
        _selectedFunction=defaultSelected;
        _minLength = 1,
        _dataField = "dataField",
        _labelField = "labelField";

    _selection=d3.select(parent);

    function component() {
    
            var container = d3.select(parent);
            
            container.attr("width", __width)
                .attr("height", __height);            
            
            var enter = container.append("div")
                    .attr("id","bp-ac")
                    .attr("class","bp-ac")
                
            var input = enter.append("input")
            			.attr("id","bp-input")
            			.attr("type","search")
                        .attr("class", "form-control")
                        .attr("placeholder",_placeHolder)
                        .attr("type","text")
                        .on("keyup",onKeyUp);

            var dropDown=enter.append("div").attr("class","bp-autocomplete-dropdown");

            var searching=dropDown.append("div").attr("class","bp-autocomplete-searching").text("Searching ...");

            hideSearching();
            hideDropDown();


            function onKeyUp() {
                _searchTerm=input.node().value;
                var e=d3.event;

                if (!(e.which == 38 || e.which == 40 || e.which == 13)) {
                    if (!_searchTerm || _searchTerm == "") {
                        showSearching("No results");
                    }
                    else if (isNewSearchNeeded(_searchTerm,_lastSearchTerm)) {
                        _lastSearchTerm=_searchTerm;
                        _currentIndex=-1;
                        _results=[];
                        showSearching();
                        search();
                        processResults();
                        if (_matches.length == 0) {
                            showSearching("No results");
                        }
                        else {
                            hideSearching();
                            showDropDown();
                        }

                    }

                }
                else {
                    e.preventDefault();
                }
            }

            function processResults() {

                var results=dropDown.selectAll(".bp-autocomplete-row").data(_matches, function (d) {
                    return d[_dataField];});
                results.enter()
                    .append("div").attr("class","bp-autocomplete-row")
                    .on("click",function (d,i) { row_onClick(d); })
                    .append("div").attr("class","bp-autocomplete-title")
                    .html(function (d) {
                        var re = new RegExp(_searchTerm, 'i');
                        var strPart = d[_dataField].match(re)[0];
                        return d[_dataField].replace(re, "<span class='bp-autocomplete-highlight'>" + strPart + "</span>");
                    });

                results.exit().remove();

                //Update results

                results.select(".bp-autocomplete-title")
                    .html(function (d,i) {
                        var re = new RegExp(_searchTerm, 'i');
                        var strPart = _matches[i][_dataField].match(re);
                        if (strPart) {
                            strPart = strPart[0];
                            return _matches[i][_dataField].replace(re, "<span class='bp-autocomplete-highlight'>" + strPart + "</span>");
                        }

                    });


            }

            function search() {

                var str=_searchTerm;
                console.log("searching on " + _searchTerm);
                console.log("-------------------");

                if (str.length >= _minLength) {
                    _matches = [];
                    for (var i = 0; i < _keys.length; i++) {
                        var match = false;
                        match = match || (_keys[i][_dataField].toLowerCase().indexOf(str.toLowerCase()) >= 0);
                        if (match && _matches.length < 10) {
                            _matches.push(_keys[i]);
                            //console.log("matches " + _keys[i][_dataField]);
                        }
                    }
                }
            }

            function row_onClick(d) {
                hideDropDown();
                input.node().value= d[_dataField];
                _selectedFunction(d);
            }

            function isNewSearchNeeded(newTerm, oldTerm) {
                return newTerm.length >= _minLength && newTerm != oldTerm;
            }

            function hideSearching() {
                searching.style("display","none");
            }

            function hideDropDown() {
                dropDown.style("display","none");
            }

            function showSearching(value) {
                searching.style("display","block");
                searching.text(value);
            }

            function showDropDown() {
                dropDown.style("display","block");
            }

    }


    function measure() {
        _width=__width - _margin.right - _margin.left;
        _height=__height - _margin.top - _margin.bottom;
    }

    function defaultSelected(d) {
        console.log(d[_dataField] + " selected");
    }


    component.render = function() {
        measure();
        component();
        return component;
    }

    component.keys = function (_) {
        if (!arguments.length) return _keys;
        _keys = _;
        return component;
    }

    component.dataField = function (_) {
        if (!arguments.length) return _dataField;
        _dataField = _;
        return component;
    }

    component.labelField = function (_) {
        if (!arguments.length) return _labelField;
        _labelField = _;
        return component;
    }

    component.margin = function(_) {
        if (!arguments.length) return _margin;
        _margin = _;
        measure();
        return component;
    };

    component.width = function(_) {
        if (!arguments.length) return __width;
        __width = _;
        measure();
        return component;
    };

    component.height = function(_) {
        if (!arguments.length) return __height;
        __height = _;
        measure();
        return component;
    };

    component.delay = function(_) {
        if (!arguments.length) return _delay;
        _delay = _;
        return component;
    };

    component.keys = function(_) {
        if (!arguments.length) return _keys;
        _keys = _;
        return component;
    };

    component.placeHolder = function(_) {
        if (!arguments.length) return _placeHolder;
        _placeHolder = _;
        return component;
    };

    component.onSelected = function(_) {
        if (!arguments.length) return _selectedFunction;
        _selectedFunction = _;
        return component;
    };



    return component;

}

var regions = ["Africa", "SouthAmerica", "WestAsia", "Oceania", "Europe", "JapanKorea", "NorthAmerica", "SoutheastAsia", "SouthAsia", "China"]

var cladeToSeq = {}

var globalDate = new Date();

var nodes, tips, rootNode, links, vaccines;

var nDisplayTips, displayRoot;
if (document.getElementById("gtspec") != null){
    var freqdefault = document.getElementById("gtspec").value;
}else{
    var freqdefault ='';
}

function treePlotHeight(width) {
	return 400 + 0.30*width;
}
var containerWidth = parseInt(d3.select(".treeplot-container").style("width"), 10);
var treeWidth = containerWidth;
var treeHeight = treePlotHeight(treeWidth);
var tree = d3.layout.tree()
	.size([treeHeight, treeWidth]);


var treeplot = d3.select("#treeplot")
	.attr("width", treeWidth)
	.attr("height", treeHeight);

var legend = d3.select("#legend")
	.attr("width", 280)
	.attr("height", 100);

var colorBy = document.getElementById("coloring").value;
var colorScale;
var epiColorDomain = genericDomain;
var nonEpiColorDomain = genericDomain;
var rbsColorDomain = genericDomain;
var dateColorDomain = genericDomain;
var dfreqColorDomain = genericDomain.map(function(d){return Math.round(100*(0.2+d*1.8))/100;});
var time_step;


d3.json(path + file_prefix + "meta.json", function(error, json) {
    if (error) return console.warn(error);
    d3.select("#updated").text(json['updated']);
    commit_id = json['commit'];
    short_id = commit_id.substring(0, 6);
    d3.select("#commit")
        .append("a")
        .attr("href", "http://github.com/blab/nextflu/commit/" + commit_id)
        .text(short_id);

});

var regions = [
["Africa",  		"#5097BA"],
["SouthAmerica",  	"#60AA9E"],
["WestAsia",  		"#75B681"],
["Oceania",  		"#8EBC66"],
["Europe",  		"#AABD52"],
["JapanKorea",  	"#C4B945"],
["NorthAmerica",  	"#D9AD3D"],
["SoutheastAsia",  	"#E59637"],
["SouthAsia",  		"#E67030"],
["China",  			"#DF4327"]
];

var ymd_format = d3.time.format("%Y-%m-%d");
var dateValues, earliestDate, dateScale, niceDateScale, counterData;

function adjust_freq_by_date() {
	calcTipCounts(rootNode);
	var tipCount = rootNode.tipCount;
	nDisplayTips = displayRoot.tipCount;
	console.log("Total tipcount: " + tipCount);
	nodes.forEach(function (d) {
		d.frequency = (d.tipCount)/tipCount;
	});
}

var drag = d3.behavior.drag()
	.on("drag", dragged)
	.on("dragstart", function() {
		d3.selectAll(".date-input-text").style("fill", "#5DA8A3");
		d3.selectAll(".date-input-marker").style("fill", "#5DA8A3");
		d3.selectAll(".date-input-window").style("stroke", "#5DA8A3");
		d3.selectAll(".date-input-edge").style("stroke", "#5DA8A3");
	})
	.on("dragend", function() {
		d3.selectAll(".date-input-text").style("fill", "#CCC");
		d3.selectAll(".date-input-marker").style("fill", "#CCC");
		d3.selectAll(".date-input-window").style("stroke", "#CCC");
		d3.selectAll(".date-input-edge").style("stroke", "#CCC");
		dragend();
	});
	
var dragMin = d3.behavior.drag()
	.on("drag", draggedMin)
	.on("dragstart", function() {
		d3.selectAll(".date-input-text").style("fill", "#5DA8A3");
		d3.selectAll(".date-input-marker").style("fill", "#5DA8A3");
		d3.selectAll(".date-input-window").style("stroke", "#5DA8A3");
		d3.selectAll(".date-input-edge").style("stroke", "#5DA8A3");
	})
	.on("dragend", function() {
		d3.selectAll(".date-input-text").style("fill", "#CCC");
		d3.selectAll(".date-input-marker").style("fill", "#CCC");
		d3.selectAll(".date-input-window").style("stroke", "#CCC");
		d3.selectAll(".date-input-edge").style("stroke", "#CCC");
		dragend();
	});	


function calcNodeAges(tw){
	tips.forEach(function (d) {
		var date = new Date(d.date);
		var oneYear = 365.25*24*60*60*1000; // days*hours*minutes*seconds*milliseconds
		var diffYears = (globalDate.getTime() - date.getTime()) / oneYear;
		d.diff = diffYears;
		if (d.diff > 0 && d.diff < tw){
			d.current  = true;
		} else{
			d.current = false;
		}
		for (var k in restrictTo){
			if (d[k]!=restrictTo[k] && restrictTo[k]!="all"){
				d.current = false;
			}
		}
	});
};



function dragged(d) {

	d.date = dateScale.invert(d3.event.x);
	d.x = dateScale(d.date);
	var startDate = new Date(d.date);
	startDate.setDate(startDate.getDate() - (time_window * 365.25));
	d.x2 = dateScale(startDate);	
	
	d3.selectAll(".date-input-text")
		.attr("dx", function(d) {return 0.5*d.x})
		.text(function(d) {
			var format = d3.time.format("%Y %b %-d");
			return format(d.date)
		});
	d3.selectAll(".date-input-marker")
		.attr("cx", function(d) {return d.x});
	d3.selectAll(".date-input-window")
		.attr("x1", function(d) {return d.x})
		.attr("x2", function(d) {return d.x2});
	d3.selectAll(".date-input-edge")
		.attr("x1", function(d) {return d.x2;})
		.attr("x2", function(d) {return d.x2});		

	globalDate = d.date;

	calcNodeAges(time_window);
	treeplot.selectAll(".link")
		.style("stroke", function(d){return "#ccc";})

	treeplot.selectAll(".tip")
		.style("visibility", tipVisibility)
		.style("fill", "#CCC")
		.style("stroke", "#AAA");

	treeplot.selectAll(".vaccine")
		.style("visibility", function(d) {
			var date = new Date(d.choice);
			var oneYear = 365.25*24*60*60*1000; // days*hours*minutes*seconds*milliseconds
			var diffYears = (globalDate.getTime() - date.getTime()) / oneYear;
			if (diffYears > 0) { return "visible"; }
				else { return "hidden"; }
			});

}

function draggedMin(d) {

	d.date = dateScale.invert(d3.event.x);
	d.x2 = dateScale(d.date);

	var oneYear = 365.25*24*60*60*1000; // days*hours*minutes*seconds*milliseconds
	time_window = (globalDate.getTime() - d.date.getTime()) / oneYear;
	
	d3.selectAll(".date-input-window")
		.attr("x2", function(d) {return d.x2});
	d3.selectAll(".date-input-edge")
		.attr("x1", function(d) {return d.x2;})
		.attr("x2", function(d) {return d.x2});		

	calcNodeAges(time_window);
	treeplot.selectAll(".link")
		.style("stroke", function(d){return "#ccc";})

	treeplot.selectAll(".tip")
		.style("visibility", tipVisibility)
		.style("fill", "#CCC")
		.style("stroke", "#AAA");

	treeplot.selectAll(".vaccine")
		.style("visibility", function(d) {
			var date = new Date(d.choice);
			var oneYear = 365.25*24*60*60*1000; // days*hours*minutes*seconds*milliseconds
			var diffYears = (globalDate.getTime() - date.getTime()) / oneYear;
			if (diffYears > 0) { return "visible"; }
				else { return "hidden"; }
			});

}

function dragend() {
	var num_date = globalDate/1000/3600/24/365.25+1970;	

	updateColorDomains(num_date);
	for (var ii=0; ii<rootNode.pivots.length-1; ii++){
		if (rootNode.pivots[ii]<num_date && rootNode.pivots[ii+1]>=num_date){
			freq_ii=Math.max(dfreq_dn,ii+1);
		}
	}
	console.log("changed frequency index to "+freq_ii+" date cut off is "+num_date);
	console.log("recalculating node ages");
	calcNodeAges(time_window);
	console.log("adjusting node colors");
	adjust_coloring_by_date();
	console.log("updating frequencies");
	adjust_freq_by_date();

	if (colorBy == "genotype") {
		colorByGenotype();
	}
	if (colorBy == "date") {
		removeLegend();
		makeLegend();
	}

	d3.selectAll(".link")
		.transition().duration(500)
		.attr("points", branchPoints)
		.style("stroke-width", branchStrokeWidth)
		.style("stroke", branchStrokeColor);				

	d3.selectAll(".tip")
		.transition().duration(500)
		.style("visibility", tipVisibility)
		.style("fill", tipFillColor)
		.style("stroke", tipStrokeColor);
				
	
	if ((typeof tip_labels != "undefined")&&(tip_labels)) {
		nDisplayTips = displayRoot.fullTipCount;
		d3.selectAll(".tipLabel")
			.transition().duration(1000)
			.style("font-size", tipLabelSize);
	}	
	
}


function date_init(){
	nodes.forEach(function (d) {d.dateval = new Date(d.date)});
	var dateValues = nodes.filter(function(d) {
		return typeof d.date === 'string';
		}).map(function(d) {
		return new Date(d.date);
	});
	earliestDate = new Date(d3.min(dateValues));
	earliestDate.setDate(earliestDate.getDate() + 1);
	
	dateScale = d3.time.scale()
		.domain([earliestDate, globalDate])
		.range([5, 205])
		.clamp([true]);	
	
	niceDateScale = d3.time.scale()
		.domain([earliestDate, globalDate])
		.range([5, 205])
		.clamp([true])
		.nice(d3.time.month);
	
	counterData = {}
	counterData['date'] = globalDate
	counterData['x'] = dateScale(globalDate)
	var startDate = new Date(globalDate);
	startDate.setDate(startDate.getDate() - (time_window * 365.25));
	counterData['x2'] = dateScale(startDate);

	d3.select("#date-input")
		.attr("width", 240)
		.attr("height", 65);

	var counter = d3.select("#date-input").selectAll(".date-input-text")
		.data([counterData])
		.enter()
		.append("text")
		.attr("class", "date-input-text")
		.attr("text-anchor", "left")
		.attr("dx", function(d) {return 0.5*d.x})		
		.attr("dy", "1.0em")
		.text(function(d) {
			var format = d3.time.format("%Y %b %-d");
			return format(d.date)
		})
		.style("cursor", "pointer")
		.call(drag);

	var customTimeFormat = d3.time.format.multi([
		[".%L", function(d) { return d.getMilliseconds(); }],
		[":%S", function(d) { return d.getSeconds(); }],
		["%I:%M", function(d) { return d.getMinutes(); }],
		["%I %p", function(d) { return d.getHours(); }],
		["%a %d", function(d) { return d.getDay() && d.getDate() != 1; }],
		["%b %d", function(d) { return d.getDate() != 1; }],
		["%b", function(d) { return d.getMonth(); }],
		["%Y", function() { return true; }]
		]);

	var dateAxis = d3.svg.axis()
		.scale(niceDateScale)
		.orient('bottom')
		.ticks(5)
		.tickFormat(customTimeFormat)
		.outerTickSize(2)
		.tickPadding(8);

	d3.select("#date-input").selectAll(".date-input-axis")
		.data([counterData])
		.enter()
		.append("g")
		.attr("class", "date-input-axis")
		.attr("transform", "translate(0,35)")
		.call(dateAxis);

	var window = d3.select("#date-input").selectAll(".date-input-window")
		.data([counterData])
		.enter()
		.append("line")
		.attr("class", "date-input-window")
		.attr("x1", function(d) { return d.x; })
		.attr("x2", function(d) { return d.x2; })	
		.attr("y1", 35)
		.attr("y2", 35)		
		.style("stroke", "#CCC")
		.style("stroke-width", 5);
		
	var edge = d3.select("#date-input").selectAll(".date-input-edge")
		.data([counterData])
		.enter()
		.append("line")
		.attr("class", "date-input-edge")
		.attr("x1", function(d) { return d.x2; })
		.attr("x2", function(d) { return d.x2; })	
		.attr("y1", 30)
		.attr("y2", 40)
		.style("stroke", "#CCC")
		.style("stroke-width", 3)	
		.style("cursor", "pointer")
		.call(dragMin);		

	var marker = d3.select("#date-input").selectAll(".date-input-marker")
		.data([counterData])
		.enter()
		.append("circle")
		.attr("class", "date-input-marker")
		.attr("cx", function(d) {return d.x})
		.attr("cy", 35)
		.attr("r", 6)
		.style("fill", "#CCC")
		.style("stroke", "#777")		
		.style("cursor", "pointer")
		.call(drag);

}
var legendRectSize = 15;
var legendSpacing = 4;
function makeLegend(){
	
	d3.select("#legend-title").text(function(d){
		if (colorBy == "ep") {
			return "Epitope mutations";
		}
		if (colorBy == "ne") {
			return "Non-epitope mutations";
		}
		if (colorBy == "rb") {
			return "Receptor binding mutations";
		}
		if (colorBy == "lbi") {
			return "Local branching index";
		}
		if (colorBy == "region") {
			return "Region";
		}
		if (colorBy == "genotype") {
			return "Genotype";
		}
		if (colorBy == "date") {
			return "Date";
		}		
		if (colorBy == "dfreq") {
			var tmp_nmonth = Math.round(12*dfreq_dn*time_step);
			var tmp_text = "Freq. change ("+tmp_nmonth+" month";
			if (tmp_nmonth>1){
				tmp_text+='s';
			}
			return tmp_text+')';
		}
	});

	// construct a dictionary that maps a legend entry to the preceding interval
	var lower_bound = {};
	lower_bound[colorScale.domain()[0]] = colorScale.domain()[0];
	for (var i=1; i<colorScale.domain().length; i++){
		lower_bound[colorScale.domain()[i]]=colorScale.domain()[i-1];
	}

	// function that equates a tip and a legend element
	// exact match is required for categorical qunantities such as genotypes, regions
	// continuous variables need to fall into the interal (lower_bound[leg], leg]
	var legend_match = function(leg, tip){
		if ((colorBy=='lbi')||(colorBy=='date')||(colorBy=='dfreq')){
			return (tip.coloring<=leg)&&(tip.coloring>lower_bound[leg]);
		}else{
			return tip.coloring==leg;
		}
	}

	var tmp_leg = legend.selectAll(".legend")
	.data(colorScale.domain())
	.enter().append('g')
	.attr('class', 'legend')
	.attr('transform', function(d, i) {
		var stack = 5;
		var height = legendRectSize + legendSpacing;
		var fromRight = Math.floor(i / stack);
		var fromTop = i % stack;
		var horz = fromRight * 145 + 5;				
		var vert = fromTop * height + 5;
		return 'translate(' + horz + ',' + vert + ')';
	 });
	tmp_leg.append('rect')
	.attr('width', legendRectSize)
	.attr('height', legendRectSize)
	.style('fill', function (d) {
	 	var col = colorScale(d);
	 	return d3.rgb(col).brighter([0.35]).toString();
	 })
	.style('stroke', function (d) {
   		var col = colorScale(d);
   		return d3.rgb(col).toString();
 	})
   .on('mouseover', function(leg){
    	treeplot.selectAll(".tip") //highlight all tips corresponding to legend
            .filter(function (d){return legend_match(leg, d);})
            .attr("r", function(d){return tipRadius*1.7;})
            .style("fill", function (t) {
              return d3.rgb(tipFillColor(t)).brighter();
            });
		}) 
  	.on('mouseout', function(leg){
    	treeplot.selectAll(".tip") //undo highlight
            .filter(function (d){return legend_match(leg, d);})
            .attr("r", function(d){return tipRadius;})
            .style("fill", function (t) {
              return d3.rgb(tipFillColor(t));
            });
	    });

	tmp_leg.append('text')
	.attr('x', legendRectSize + legendSpacing + 5)
	.attr('y', legendRectSize - legendSpacing)
	.text(function(d) {
		var label = d.toString().replace(/([a-z])([A-Z])/g, '$1 $2').replace(/,/g, ', ');
		if (colorBy == "dfreq") {
			label += "\u00D7";
		}
		return label;
	})
   .on('mouseover', function(leg){
    	treeplot.selectAll(".tip")
            .filter(function (d){return legend_match(leg, d);})
            .attr("r", function(d){return tipRadius*1.7;})
            .style("fill", function (t) {
              return d3.rgb(tipFillColor(t)).brighter();
            });
		}) 
  	.on('mouseout', function(leg){
    	treeplot.selectAll(".tip")
            .filter(function (d){return legend_match(leg, d);})
            .attr("r", function(d){return tipRadius;})
            .style("fill", function (t) {
              return d3.rgb(tipFillColor(t));
            });
	    });
	return tmp_leg;
}

function removeLegend(){
	legend.selectAll('.legend')
  .remove();
}


// 2 color	["#5097BA", "#DF4327"]

var colors = [
	[],
	["#8EBC66"],
	["#5EA9A1", "#E0A23A"],
	["#5CA7A4", "#A4BE56", "#E39B39"],
	["#5AA5A7", "#8ABB6A", "#C3BA46", "#E69237"],
	["#59A3AA", "#7CB879", "#AABD52", "#D4B13F", "#E68634"],
	["#57A1AD", "#73B583", "#98BD5E", "#C0BA47", "#DEA73C", "#E67B32"],
	["#559EB1", "#6DB28B", "#8CBB68", "#AFBD4F", "#CFB541", "#E39B39", "#E56E30"],
	["#539CB4", "#68AF92", "#83BA71", "#A2BE57", "#C1BA47", "#DAAC3D", "#E68E36", "#E35F2D"],
	["#519AB7", "#64AC99", "#7BB879", "#97BD5F", "#B5BD4C", "#CFB541", "#E0A13A", "#E67F33", "#E1512A"],
	["#5097BA", "#60AA9E", "#75B681", "#8EBC66", "#AABD52", "#C4B945", "#D9AD3D", "#E59637", "#E67030", "#DF4327"],
	["#4E95BD", "#5DA8A3", "#70B487", "#87BB6D", "#A0BE58", "#BABC4A", "#D0B440", "#E0A23A", "#E68634", "#E35F2D", "#DD3425"],
	["#4D92BF", "#5AA5A8", "#6BB18D", "#80B974", "#98BD5E", "#B1BD4E", "#C8B944", "#DAAC3D", "#E59738", "#E67732", "#E14F2A", "#DB2522"]
];
var genotypeColors = ["#60AA9E", "#D9AD3D", "#5097BA", "#E67030", "#8EBC66", "#E59637", "#AABD52", "#DF4327", "#C4B945", "#75B681"];

var epitopeColorScale = d3.scale.linear().clamp([true])
	.domain(epiColorDomain)
	.range(colors[10]);

var nonepitopeColorScale = d3.scale.linear().clamp([true])
	.domain(nonEpiColorDomain)
	.range(colors[10]);

var receptorBindingColorScale = d3.scale.linear().clamp([true])
	.domain(rbsColorDomain)
	.range(colors[4]);

var lbiColorScale = d3.scale.linear()
	.domain([0.0, 0.02, 0.04, 0.07, 0.1, 0.2, 0.4, 0.7, 0.9, 1.0])
	.range(colors[10]);

var dfreqColorScale = d3.scale.linear()
	.domain(dfreqColorDomain)
	.range(colors[10]);

var dateColorScale = d3.scale.linear().clamp([true])
	.domain(dateColorDomain)
	.range(colors[10]);

var regionColorScale = d3.scale.ordinal()
	.domain(regions.map(function(d){return d[0];}))
	.range(regions.map(function(d){return d[1];}));

// "ep", "ne" and "rb" need no adjustments
function adjust_coloring_by_date() {
	if (colorBy == "lbi") {
		calcLBI(rootNode, nodes, false);
		nodes.forEach(function (d) {
			d.coloring = d.LBI;
		});
	}
	else if (colorBy == "dfreq") {
		calcDfreq(rootNode, freq_ii);
		nodes.forEach(function (d) {
			d.coloring = d.dfreq;
		});
	}
	else if (colorBy == "date") {
		nodes.forEach(function (d) {
			d.coloring = d.num_date;
		});
	}
}

function stateAtPosition(clade, gene, pos){
	if (typeof cladeToSeq[clade][gene][pos] == "undefined"){
		return cladeToSeq["root"][gene][pos];
	}else{
		return cladeToSeq[clade][gene][pos];
	}
}

function colorByTrait() {

	colorBy = document.getElementById("coloring").value;
	console.log(colorBy);

	if (colorBy == "ep") {
		colorScale = epitopeColorScale;
		nodes.map(function(d) { d.coloring = d.ep; });
	}
	else if (colorBy == "ne") {
		colorScale = nonepitopeColorScale;
		nodes.map(function(d) { d.coloring = d.ne; });
	}
	else if (colorBy == "rb") {
		colorScale = receptorBindingColorScale;
		nodes.map(function(d) { d.coloring = d.rb; });
	}
	else if (colorBy == "lbi") {
		colorScale = lbiColorScale;
		adjust_coloring_by_date();
	}
	else if (colorBy == "dfreq") {
		colorScale = dfreqColorScale;
		adjust_coloring_by_date();
	}
	else if (colorBy == "region") {
		colorScale = regionColorScale;
		nodes.map(function(d) { d.coloring = d.region; });
	}
	else if (colorBy == "date") {
		colorScale = dateColorScale;
		nodes.map(function(d) { d.coloring = d.num_date; });
	}

	treeplot.selectAll(".link")
		.style("stroke", branchStrokeColor);

	d3.selectAll(".tip")
		.style("visibility", tipVisibility)
		.style("fill", tipFillColor)
		.style("stroke", tipStrokeColor);

	if (typeof tree_legend != undefined){
		removeLegend();
	}
	tree_legend = makeLegend();
}

function tipStrokeColor(d) {
	var col = colorScale(d.coloring);
	return d3.rgb(col).toString();
}

function tipFillColor(d) {
	var col = colorScale(d.coloring);	;
	return d3.rgb(col).brighter([0.65]).toString();
}

function branchStrokeColor(d) {
	var col;
	if (colorBy == "region" || colorBy == "date") {
		col = "#AAA";
	}
	else {
		col = colorScale(d.target.coloring);
	}
	var modCol = d3.interpolateRgb(col, "#BBB")(0.6);
	return d3.rgb(modCol).toString();
}

function colorByGenotype() {
	var positions_string = document.getElementById("gt-color").value.split(',');
	var positions_list = []
	positions_string.map(function(d) {
		var pos_fields = d.split(':');
		var val, gene;
		if (pos_fields.length==1){
			val = parseInt(pos_fields[0])-1;
			gene=default_gene;
		}else if (pos_fields.length==2){
			val = parseInt(pos_fields[1])-1;
			gene=pos_fields[0].replace(' ','');
		}else{
			val = parseInt('NaN');
		}
		console.log('attempt genotype coloring: '+ [gene, val]);
		if ((!isNaN(val))&&(typeof cladeToSeq["root"][gene]!="undefined")) {
			if (val < cladeToSeq["root"][gene].length) {
				positions_list.push([gene, val]);
			}
		}
	});
	console.log(positions_list);
	if (positions_list.length > 0) {
		colorBy = "genotype";
		colorByGenotypePosition(positions_list);
	}
	else {
		d3.select("#coloring").each(colorByTrait);
		gt = parse_gt_string(freqdefault);
		make_gt_chart(gt);
		document.getElementById("gtspec").value = freqdefault;
	}
}

function colorByGenotypePosition (positions) {
	var gts = nodes.map(function (d) {
		var tmp = [];
		for (var i=0; i<positions.length; i++){
			tmp[tmp.length] = positions[i][0]+':'+(positions[i][1]+1)+stateAtPosition(d.clade, positions[i][0], positions[i][1]);
		}
		d.coloring = tmp.join('/');
		return d.coloring;});
	var unique_gts = d3.set(gts).values();
	var gt_counts = {};
	for (var i=0; i<unique_gts.length; i++){gt_counts[unique_gts[i]]=0;}
	gts.forEach(function (d) {gt_counts[d]+=1;});
	var filtered_gts = unique_gts.filter(function (d) {return gt_counts[d]>=10;});
	filtered_gts.sort(function (a,b){
		var res;
		if (gt_counts[a]>gt_counts[b]){ res=-1;}
		else if (gt_counts[a]<gt_counts[b]){ res=1;}
		else {res=0;}
		return res;});
	console.log("genotypes passed filtering:"+filtered_gts);
	colorScale = d3.scale.ordinal()
		.domain(filtered_gts)
		.range(genotypeColors);
	treeplot.selectAll(".link")
		.style("stroke", branchStrokeColor);
	treeplot.selectAll(".tip")
		.style("fill", tipFillColor)
		.style("stroke", tipStrokeColor);
	if (typeof tree_legend != undefined){
		removeLegend();
	}
	tree_legend = makeLegend();

	if ((positions.length==1)&&(filtered_gts.length>1)){
		var tmp_gts=[];
		for (var ii=0; ii<filtered_gts.length; ii+=1){
			tmp_gts.push(["global", filtered_gts[ii]])
		}
		make_gt_chart(tmp_gts);
		document.getElementById("gtspec").value = tmp_gts.map( function (d) {return d[1];}).join(', ');
	}
}

d3.select("#coloring")
	.style("cursor", "pointer")
	.on("change", colorByTrait);


var genotypeColoringEvent;
d3.select("#gt-color")
	.on("keyup", function(){
		if (typeof genotypeColoringEvent != "undefined"){clearTimeout(genotypeColoringEvent);}
		genotypeColoringEvent = setTimeout(colorByGenotype, 200);
	});






function gatherTips(node, tips) {
	if (typeof node.children != "undefined") {
		for (var i=0, c=node.children.length; i<c; i++) {
			gatherTips(node.children[i], tips);
		}
	}
	else {
		tips.push(node);
	}
	return tips;
}

function getVaccines(tips) {
	vaccines = [];
	tips.forEach(function (tip) {
		if (vaccineStrains.indexOf(tip.strain) != -1) {
			tip.choice = vaccineChoice[tip.strain];
			vaccines.push(tip);
		}
	})
	return vaccines;
}

function minimumAttribute(node, attr, min) {
	if (typeof node.children != "undefined") {
		for (var i=0, c=node.children.length; i<c; i++) {
			min = minimumAttribute(node.children[i], attr, min);
		}
	}
	else {
		if (node[attr] < min) {
			min = node[attr];
		}
	}
	return min;
}

function maximumAttribute(node, attr, max) {
	if (typeof node.children != "undefined") {
		for (var i=0, c=node.children.length; i<c; i++) {
			max = maximumAttribute(node.children[i], attr, max);
		}
	}
	else {
		if (node[attr] > max) {
			max = node[attr];
		}
	}
	return max;
}

function calcBranchLength(node){
	if (typeof node.children != "undefined") {
	for (var i=0, c=node.children.length; i<c; i++) {
		calcBranchLength(node.children[i]);
		node.children[i].branch_length = node.children[i].xvalue-node.xvalue;
	}
	}
};

/**
 * for each node, calculate the number of subtending tips (alive or dead) 
**/
function calcFullTipCounts(node){
	node.fullTipCount = 0;
	if (typeof node.children != "undefined") {
		for (var i=0; i<node.children.length; i++) {
			calcFullTipCounts(node.children[i]);
			node.fullTipCount += node.children[i].fullTipCount;
		}
	}
	else { 
		node.fullTipCount = 1;
	}
};

/**
 * for each node, calculate the number of tips in the currently selected time window. 
**/
function calcTipCounts(node){
	node.tipCount = 0;
	if (typeof node.children != "undefined") {
		for (var i=0; i<node.children.length; i++) {
			calcTipCounts(node.children[i]);
			node.tipCount += node.children[i].tipCount;
		}
	}
	else if (node.current){ 
		node.tipCount = 1;
	}
};

/**
sets each node in the tree to alive=true if it has at least one descendent with current=true
**/
function setNodeAlive(node){
	if (typeof node.children != "undefined") {
		var aliveChildren=false;
		for (var i=0, c=node.children.length; i<c; i++) {
			setNodeAlive(node.children[i]);
			aliveChildren = aliveChildren||node.children[i].alive
		}   
		node.alive = aliveChildren;
	}else{
		node.alive = node.current;
	}
};

/**
 * for each node, calculate the exponentially attenuated tree length below the node
 * the polarizer is send "up", i.e. to parents
**/
function calcUpPolarizers(node){
	node.up_polarizer = 0;
	if (typeof node.children != "undefined") {
		for (var i=0; i<node.children.length; i++) {
		calcUpPolarizers(node.children[i]);
		node.up_polarizer += node.children[i].up_polarizer;
		}
	}
	bl =  node.branch_length/LBItau;
	node.up_polarizer *= Math.exp(-bl);
	if (node.alive){ // only alive branches contribute anything
		node.up_polarizer += LBItau*(1-Math.exp(-bl));
	}
};

/**
 * for each node, calculate the exponentially attenuated tree length above the node,
 * that is "outside" the clade defined by this node. this down polarizer is send to children
**/
function calcDownPolarizers(node){
	if (typeof node.children != "undefined") {
	for (var i1=0; i1<node.children.length; i1++) {
		node.children[i1].down_polarizer = node.down_polarizer;
		for (var i2=0; i2<node.children.length; i2++) {
			if (i1!=i2){
			node.children[i1].down_polarizer += node.children[i2].up_polarizer;
			}
		}
		// account for the attenuation over the branch_length 
		bl =  node.children[i1].branch_length/LBItau;
		node.children[i1].down_polarizer *= Math.exp(-bl);
		if (node.children[i1].alive) { //the branch contributes only when the node is alive 
			node.children[i1].down_polarizer += LBItau*(1-Math.exp(-bl));
		}
		calcDownPolarizers(node.children[i1]);
	}
	}
};

function calcPolarizers(node){
	calcUpPolarizers(node);
	node.down_polarizer = 0; // set the down polarizer of the root to 0
	calcDownPolarizers(node);
};

/**
 * calculate the LBI for all nodes downstream of node
 * allnodes is provided for easy normalization at the end
**/
function calcLBI(node, allnodes){
	setNodeAlive(node);
	calcPolarizers(node);
	allnodes.forEach(function (d) {
		d.LBI=0;
		d.LBI+=d.down_polarizer;
		if (typeof d.children != "undefined") {
			for (var i=0; i<d.children.length; i++) {
				d.LBI += d.children[i].up_polarizer;
			}
		}
	});
	// normalize the LBI to range [0,1]
	maxLBI = d3.max(allnodes.map(function (d) {return d.LBI;}));
	allnodes.forEach(function (d){ d.LBI /= maxLBI;});
};

var virusTooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0, 12])
	.html(function(d) {
	
		string = "";
				
		// safe to assume the following attributes
		if (typeof d.strain != "undefined") {
			string += d.strain;
		}
		string += "<div class=\"smallspacer\"></div>";
		
		string += "<div class=\"smallnote\">";		
		
		// check if vaccine strain
		if (vaccineStrains.indexOf(d.strain) != -1) {
			string += "Vaccine strain<br>";
			var vaccine_date = new Date(vaccineChoice[d.strain]);

			string += "First chosen " + vaccine_date.toLocaleString("en-us", { month: "short" }) + " " + vaccine_date.getFullYear() + "<br>";
			string += "<div class=\"smallspacer\"></div>";
		}			
		
		if (typeof d.country != "undefined") {
			string += d.country.replace(/([A-Z])/g, ' $1');
		}
		if (typeof d.date != "undefined") {
			string += ", " + d.date;
		}
		if ((typeof d.db != "undefined") && (typeof d.accession != "undefined") && (d.db == "GISAID")) {
			string += "<br>GISAID ID: EPI" + d.accession;
		}
		if (typeof d.lab != "undefined") {
			if (d.lab != "") {
				string += "<br>Source: " + d.lab.substring(0,25);
				if (d.lab.length>25) string += '...';
			}
		}			
		string += "</div>";
		
		string += "<div class=\"smallspacer\"></div>";
				
		// following may or may not be present
		string += "<div class=\"smallnote\">";
		if (typeof d.ep != "undefined") {
			string += "Epitope distance: " + d.ep + "<br>";
		}
		if (typeof d.ne != "undefined") {
			string += "Non-epitope distance: " + d.ne + "<br>";
		}
		if (typeof d.rb != "undefined") {
			string += "Receptor binding distance: " + d.rb + "<br>";
		}
		if (typeof d.LBI != "undefined") {
			string += "Local branching index: " + d.LBI.toFixed(3) + "<br>";
		}
		string += "</div>";
		return string;
	});
treeplot.call(virusTooltip);


var linkTooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0, 12])
	.html(function(d) {
		string = ""
		if (typeof d.frequency != "undefined") {
			string += "Frequency: " + (100 * d.frequency).toFixed(1) + "%"
		}
		string += "<div class=\"smallspacer\"></div>";
		string += "<div class=\"smallnote\">";
		if ((typeof d.aa_muts !="undefined")&&(mutType=='aa')){
			var ncount = 0;
			for (tmp_gene in d.aa_muts) {ncount+=d.aa_muts[tmp_gene].length;}
			if (ncount) {string += "<b>Mutations:</b><ul>";}
			for (tmp_gene in d.aa_muts){
				if (d.aa_muts[tmp_gene].length){
					string+="<li>"+tmp_gene+":</b> "+d.aa_muts[tmp_gene].replace(/,/g, ', ') + "</li>";
				}
			}
		}
		else if ((typeof d.nuc_muts !="undefined")&&(mutType=='nuc')&&(d.nuc_muts.length)){
			var tmp_muts = d.nuc_muts.split(',');
			var nmuts = tmp_muts.length;
			tmp_muts = tmp_muts.slice(0,Math.min(10, nmuts))
			string += "<li>"+tmp_muts.join(', ');
			if (nmuts>10) {string+=' + '+ (nmuts-10) + ' more';}
			string += "</li>";
		}
		string += "</ul>";
		string += "click to zoom into clade"
		string += "</div>";
		return string;
	});
treeplot.call(linkTooltip);


var matchTooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0, 12])
	.html(function(d) {
		string = d.strain+ "<i> is closest match of:</i><ul>";
		string += "<div class=\"smallspacer\"></div>";
		for (var mi=0; mi<d.matches.length;mi++){
			string+="<li>" +d.matches[mi].substring(0,Math.min(30,d.matches[mi].length))+'</li>';
		}
		string += "</ul>";
		return string;
	});
treeplot.call(matchTooltip);

console.log('Enter tree.js');

var freqScale = d3.scale.sqrt()
	.domain([0, 1])
	.range([1, 10]);

var tipRadius = 4.0;
var left_margin = 10;
var bottom_margin = 10;
var top_margin = 32;
if ((typeof branch_labels != "undefined")&&(branch_labels)) {top_margin +=15;}
var right_margin = 10;

function initDateColorDomain(intAttributes){
	var numDateValues = tips.map(function(d) {return d.num_date;})
	var minDate = d3.min(numDateValues.filter(function (d){return d!="undefined";}));
	var maxDate = d3.max(numDateValues.filter(function (d){return d!="undefined";}));	
	if (typeof time_window == "undefined"){
		time_window = maxDate-minDate;
		console.log("defining time window as " + time_window);
	} 
	if (time_window>1){
		dateColorDomain = genericDomain.map(function (d){return Math.round(10*(maxDate - (1.0-d)*time_window))/10;});
	}else{
		dateColorDomain = genericDomain.map(function (d){return Math.round(100*(maxDate - (1.0-d)*time_window))/100;});		
	}
	console.log('setting date domain '+dateColorDomain);
	dateColorScale.domain(dateColorDomain);
}


function initColorDomain(attr, tmpCS){
	//var vals = tips.filter(function(d) {return tipVisibility(d)=='visible';}).map(function(d) {return d[attr];});
	var vals = tips.map(function(d) {return d[attr];});
	var minval = d3.min(vals);
	var maxval = d3.max(vals);	
	var rangeIndex = Math.min(10, maxval - minval + 1);
	var domain = [];
	if (maxval-minval<20)
	{
		for (var i=maxval - rangeIndex + 1; i<=maxval; i+=1){domain.push(i);}
	}else{
		for (var i=1.0*minval; i<=maxval; i+=(maxval-minval)/9.0){domain.push(i);}		
	}
	tmpCS.range(colors[rangeIndex]);
	tmpCS.domain(domain);
}

function updateColorDomains(num_date){
	dateColorDomain = genericDomain.map(function(d) {return Math.round(10*(num_date - time_window*(1.0-d)))/10;});
	dateColorScale.domain(dateColorDomain);
}

function tipVisibility(d) {
	if ((d.diff < 0 || d.diff > time_window)&(date_select==true)) {
		return "hidden";
	}
	for (var k in restrictTo){
		if (d[k]!=restrictTo[k] && restrictTo[k]!="all"){
			return "hidden";
		}
	}
	return "visible";
}

function branchPoints(d) {
	var mod = 0.5 * freqScale(d.target.frequency) - freqScale(0);
	return (d.source.x-mod).toString() + "," + d.source.y.toString() + " "
		+ (d.source.x-mod).toString() + "," + d.target.y.toString() + " "
		+ (d.target.x).toString() + "," + d.target.y.toString();
}

function branchStrokeWidth(d) {
	return freqScale(d.target.frequency);
}

function branchLabelText(d) {
	var tmp_str = d.aa_muts.replace(/,/g, ', '); 
	if (tmp_str.length>50){
		return tmp_str.substring(0,45)+'...';
	}
	else {
		return tmp_str;
	}
}

function tipLabelText(d) {
	if (d.strain.length>32){
		return d.strain.substring(0,30)+'...';
	}
	else {
		return d.strain;
	}
}

function branchLabelSize(d) {
	var n = nDisplayTips;
	if (d.fullTipCount>n/15) {
		return "12px";
	}
	else {
		return "0px";
	}
}

function tipLabelSize(d) {
	if (tipVisibility(d)!="visible"){
		return 0;
	}
	var n = nDisplayTips;
	if (n<25){
		return 16;
	}else if (n<50){
		return 12;
	}else if (n<75){
		return 8;
	}
	else {
		return 0;
	}
}

function tipLabelWidth(d) {
	return tipLabelText(d).length * tipLabelSize(d) * 0.5;
}

function tree_init(){
	calcFullTipCounts(rootNode);
	calcBranchLength(rootNode);
	rootNode.branch_length= 0.01;	
	rootNode.dfreq = 0.0;
	if (typeof rootNode.pivots != "undefined"){
		time_step = rootNode.pivots[1]-rootNode.pivots[0];		
	}else{
		time_step = 1.0/12;
	}
	//setting index of frequency trajectory to use for calculating frequency change
	freq_ii = 1;
	if (typeof rootNode.pivots != "undefined") {
		if (typeof rootNode.pivots.length != "undefined") {
			freq_ii = rootNode.pivots.length - 1;
		}
	}
	calcNodeAges(LBItime_window);
	colorByTrait();
	adjust_freq_by_date();
	tree_legend = makeLegend();
	nDisplayTips = displayRoot.fullTipCount;	
}

d3.json(path + file_prefix + "tree.json", function(error, root) {

	if (error) return console.warn(error);

	nodes = tree.nodes(root);
	links = tree.links(nodes);
	var tree_legend;
	rootNode = nodes[0];
	displayRoot = rootNode;
	tips = gatherTips(rootNode, []);
	vaccines = getVaccines(tips);

	initDateColorDomain();
	if (typeof rootNode['ep'] != "undefined"){ initColorDomain('ep', epitopeColorScale);}
	if (typeof rootNode['ne'] != "undefined"){ initColorDomain('ne', nonepitopeColorScale);}
	if (typeof rootNode['rb'] != "undefined"){ initColorDomain('rb', receptorBindingColorScale);}
	date_init();
	tree_init();

	var xValues = nodes.map(function(d) {
		return +d.xvalue;
	});

	var yValues = nodes.map(function(d) {
		return +d.yvalue;
	});

	var clade_freq_event;
	var link = treeplot.selectAll(".link")
		.data(links)
		.enter().append("polyline")
		.attr("class", "link")
		.style("stroke-width", branchStrokeWidth)
		.style("stroke", branchStrokeColor)		
		.style("stroke-linejoin", "round")
		.style("cursor", "pointer")
		.on('mouseover', function (d){
			linkTooltip.show(d.target, this);
			if ((colorBy!="genotype")&(typeof addClade !="undefined")){
				clade_freq_event = setTimeout(addClade, 1000, d);
			}
			})
		.on('mouseout', function(d) {
			linkTooltip.hide(d);
			if (typeof addClade !="undefined") {clearTimeout(clade_freq_event);};})		
		.on('click', zoom);

	if ((typeof branch_labels != "undefined")&&(branch_labels)){
		var mutations = treeplot.selectAll(".branchLabel")
			.data(nodes)
			.enter()
			.append("text")
			.attr("class", "branchLabel")
			.style("font-size", branchLabelSize)			
			.style("text-anchor", "end")
			.text(branchLabelText);
		}

	if ((typeof tip_labels != "undefined")&&(tip_labels)){
		treeplot.selectAll(".tipLabel").data(tips)
			.enter()
			.append("text")
			.attr("class","tipLabel")
			.style("font-size", function(d) {console.log(tipLabelText(d)); return tipLabelSize(d)+"px"; })
			.text(tipLabelText);
	}

	var tipCircles = treeplot.selectAll(".tip")
		.data(tips)
		.enter()
		.append("circle")
		.attr("class", "tip")
		.attr("id", function(d) { return (d.strain).replace(/\//g, ""); })
		.attr("r", tipRadius)
		.style("visibility", tipVisibility)
		.style("fill", tipFillColor)
		.style("stroke", tipStrokeColor)
		.on('mouseover', function(d) {
			virusTooltip.show(d, this);
		})
		.on('click', function(d) {
			if ((typeof d.db != "undefined") && (d.db == "GISAID") && (typeof d.accession != "undefined")) {
				var url = "http://gisaid.org/EPI/"+d.accession;
				console.log("opening url "+url);
				var win = window.open(url, '_blank');
  				win.focus();
  			}	
  		})		
		.on('mouseout', virusTooltip.hide);

	var vaccineCircles = treeplot.selectAll(".vaccine")
		.data(vaccines)
		.enter()
		.append("text")
		.attr("class", "vaccine")
		.attr('text-anchor', 'middle')
		.attr('dominant-baseline', 'central')
		.style("font-size", "28px")
		.style('font-family', 'FontAwesome')
		.style("fill", "#555555")
		.text(function(d) { return '\uf00d'; })
		.style("cursor", "default")
		.on('mouseover', function(d) {
			virusTooltip.show(d, this);
		})
		.on('mouseout', virusTooltip.hide);


	d3.select("#reset")
		.on("click", function(d) {
			displayRoot = rootNode;
			nDisplayTips = displayRoot.fullTipCount;
			var dMin = d3.min(xValues),
				dMax = d3.max(xValues),
				lMin = d3.min(yValues),
				lMax = d3.max(yValues);
			rescale(dMin, dMax, lMin, lMax);
			removeClade();
		})

	/*
	 * zoom into the tree upon click onto a branch
	 */
	function zoom(d){
		if ((colorBy!="genotype")&(typeof addClade !="undefined")){
			addClade(d);
		}
		var dy = yScale.domain()[1]-yScale.domain()[0];
		displayRoot = d.target;
		var dMin = 0.5 * (minimumAttribute(d.target, "xvalue", d.target.xvalue) + minimumAttribute(d.source, "xvalue", d.source.xvalue)),
			dMax = maximumAttribute(d.target, "xvalue", d.target.xvalue),
			lMin = minimumAttribute(d.target, "yvalue", d.target.yvalue),
			lMax = maximumAttribute(d.target, "yvalue", d.target.yvalue);
		if (dMax == dMin || lMax == lMin) {
			displayRoot = d.source;
			dMin = minimumAttribute(d.source, "xvalue", d.source.xvalue),
			dMax = maximumAttribute(d.source, "xvalue", d.source.xvalue),
			lMin = minimumAttribute(d.source, "yvalue", d.source.yvalue),
			lMax = maximumAttribute(d.source, "yvalue", d.source.yvalue);			
		}
		if ((lMax-lMin)>0.999*dy){
			lMin = lMax - dy*0.7 
		}
		var visibleXvals = tips.filter(function (d){return (d.yvalue>=lMin)&&(d.yvalue<lMax)}).map(function(d){return +d.xvalue;});
		nDisplayTips = visibleXvals.length;
		dMax = Math.max.apply(Math, visibleXvals);
		console.log("nodes in view: "+nDisplayTips+' max Xval: '+dMax);
		rescale(dMin, dMax, lMin, lMax);
	}

	/*
	 * adjust margins such that the tip labels show
	 */
	function setMargins(){
		containerWidth = parseInt(d3.select(".treeplot-container").style("width"), 10);
		treeWidth = containerWidth;
		treeHeight = treePlotHeight(treeWidth);			
		d3.select("#treeplot")
			.attr("width", treeWidth)
			.attr("height", treeHeight);
		if ((typeof tip_labels != "undefined")&&(tip_labels)){
			var maxTextWidth = 0;
			var labels = treeplot.selectAll(".tipLabel")
				.data(tips)
				.each(function(d) {
					var textWidth = tipLabelWidth(d);
					if (textWidth>maxTextWidth) {
						maxTextWidth = textWidth;
					}
				});
			right_margin = maxTextWidth + 10;
		}							
		xScale.range([left_margin, treeWidth - right_margin]);
		yScale.range([top_margin, treeHeight - bottom_margin]);
	}

	/*
	 * rescale the tree to a window defined by the arguments
	 * dMin, dMax  -- minimal and maximal horizontal dimensions
	 * lMin, lMax  -- minimal and maximal vertical dimensions
	 */
	function rescale(dMin, dMax, lMin, lMax) {
		setMargins();
		xScale.domain([dMin,dMax]);
		yScale.domain([lMin,lMax]);
		virusTooltip.hide();
		linkTooltip.hide();
		matchTooltip.hide();				
		transform(1500)
	}

	/*
	 *move all svg items to their new location
	 *dt -- the duration of the transition
	 */
	function transform(dt){
		nodes.forEach(function (d) {
			d.x = xScale(d.xvalue);
			d.y = yScale(d.yvalue);
		});

		treeplot.selectAll(".tip")
			.transition().duration(dt)
			.attr("cx", function(d) { return d.x; })
			.attr("cy", function(d) { return d.y; });

		treeplot.selectAll(".vaccine")
			.transition().duration(dt)
			.attr("x", function(d) { return d.x; })
			.attr("y", function(d) { return d.y; });

		treeplot.selectAll(".seqmatch")
			.transition().duration(dt)
			.attr("x", function(d) { return d.x; })
			.attr("y", function(d) { return d.y; });
			
		treeplot.selectAll(".strainmatch")
			.transition().duration(dt)
			.attr("x", function(d) { return d.x; })
			.attr("y", function(d) { return d.y; });			

		treeplot.selectAll(".link")
			.transition().duration(dt)
			.attr("points", branchPoints);

		if ((typeof tip_labels != "undefined")&&(tip_labels)){
			treeplot.selectAll(".tipLabel")
				.transition().duration(dt)
				.style("font-size", function(d) {return tipLabelSize(d)+"px"; })			
				.attr("x", function(d) { return d.x+10; })
				.attr("y", function(d) { return d.y+4; });
		}	
			
		if ((typeof branch_labels != "undefined")&&(branch_labels)){
			console.log('shift branch_labels');
			treeplot.selectAll(".branchLabel")
				.transition().duration(dt)
				.style("font-size", branchLabelSize)				
				.attr("x", function(d) {  return d.x - 6;})
				.attr("y", function(d) {  return d.y - 3;});
		}

		if (typeof clades !="undefined"){
			treeplot.selectAll(".annotation")
				.transition().duration(dt)
				.attr("x", function(d) {
					return xScale(d[1]) - 10;
				})
				.attr("y", function(d) {
					return yScale(d[2]) - 6;
				});
		}
	}

	d3.select(window).on('resize', resize); 
	
	function resize() {
		setMargins();
		transform(0);
	}
	

	function restrictToFunc(rt) {
		restrictTo[rt] = document.getElementById(rt).value;
		console.log("restriction to "+rt+" "+restrictTo[rt]);	
		d3.selectAll(".tip")
			.style("visibility", tipVisibility);
		dragend();		
	}

	for (rt in restrictTo){
		var tmp = document.getElementById(rt);
		if (tmp!=null){
			restrictTo[rt] = tmp.value;
		}else{restrictTo[rt]='all';}
		console.log(restrictTo);
		d3.select("#"+rt)
			.style("cursor", "pointer")
			.on("change", (function(restrictor){
							return function(){
								return restrictToFunc(restrictor);
							}
						})(rt));		
	}


	var mc = autocomplete(document.getElementById('straininput'))
		.keys(tips)
		.dataField("strain")
		.placeHolder("search strains...")
		.width(800)
		.height(500)
		.onSelected(highlightStrainSearch)
		.render();


	// add clade labels
	clades = rootNode["clade_annotations"];
	if (typeof clades != "undefined"){
		console.log(clades);
		var clade_annotations = treeplot.selectAll('.annotation')
			.data(clades)
			.enter()
			.append("text")
			.attr("class", "annotation")
			.style("text-anchor", "end")
			.text(function (d) {
				return d[0];
			});
		}
	var xScale = d3.scale.linear()
		.domain([d3.min(xValues), d3.max(xValues)]);
	var yScale = d3.scale.linear()
		.domain([d3.min(yValues), d3.max(yValues)]);
	setMargins();

	resize();
});

d3.json(path + file_prefix + "sequences.json", function(error, json) {
	if (error) return console.warn(error);
	cladeToSeq=json;
});


var frequencies, pivots;
var gene = 'nuc';
var mutType = 'aa';
/**
 * for each node, calculate the derivative of the frequency tranjectory. if none exists, copy parent
**/
function calcDfreq(node, freq_ii){
	if (typeof node.children != "undefined") {
		for (var i1=0; i1<node.children.length; i1++) {
			if (typeof node.children[i1].freq != "undefined") {
				if (node.children[i1].freq["global"] != "undefined"){
					var tmp_freq = node.children[i1].freq["global"]
					//node.children[i1].dfreq = 0.5*(tmp_freq[freq_ii] - tmp_freq[freq_ii-dfreq_dn])/(tmp_freq[freq_ii] + tmp_freq[freq_ii-dfreq_dn] + 0.2);
					node.children[i1].dfreq = (tmp_freq[freq_ii] + 0.01)/(tmp_freq[freq_ii-dfreq_dn] + 0.03);
				} else {
					node.children[i1].dfreq = node.dfreq;
				}
			}
			calcDfreq(node.children[i1], freq_ii);
		}
	}
};
function parse_gt_string(gt){
	mutations = [];
	gt.split(',').map( function (d) {
		var tmp = d.split(/[\s//]/); //FIXME: make more inclusive
		var region;
		var positions = [];
		for (var i=0; i<tmp.length; i++){
			if (contains(["EU","NA","AS","OC"], tmp[i])){
				region = tmp[i];
			}else{
				if (tmp[i].length>0) positions.push(tmp[i]);
			}
		}
		if (typeof region == "undefined") region="global";
		// sort if this is a multi mutation genotype
		if (positions.length>1){
			positions.sort(function (a,b){
				return parseInt(a.substring(0,a.length-1)) - parseInt(b.substring(0,b.length-1));
			});
		}
		mutations.push([region, positions.join('/')]);
	});
	return mutations;
};

/**
loops over all genotypes from a certain region and sums the frequency contributions
of the genotype matches at the specified positions
**/
function get_frequencies(region, gt){
	var freq = [];
	for (var pi=0; pi<pivots.length; pi++){freq[freq.length]=0;}
	console.log("searching for "+region+' ' + gt);
	if (frequencies["clades"][region][gt.toLowerCase()]!=undefined) {
		console.log(gt+" found as clade");
		for (var pi=0; pi<freq.length; pi++){
			freq[pi]+=frequencies["clades"][region][gt.toLowerCase()][pi];
		}
	}else if ((typeof frequencies["genotypes"] !="undefined") && (frequencies["genotypes"][region][gt]!=undefined)) {
		console.log(gt+" found as genotype");
		for (var pi=0; pi<freq.length; pi++){
			freq[pi]+=frequencies["genotypes"][region][gt][pi];
		}
	}else if (frequencies["mutations"] !== undefined){
		var tmp_mut = gt.split(':');
		var mut ="";
		if (tmp_mut.length==1){
			mut = default_gene+":"+gt;
		}else{
			mut = gt;
		}
		if (frequencies["mutations"][region][mut]!==undefined) {
			console.log(gt+" found as mutation");
			for (var pi=0; pi<freq.length; pi++){
				freq[pi]+=frequencies["mutations"][region][mut][pi];
			}
		}else{
			console.log("not found "+gt);
		}
	}else{
		console.log("not found "+gt);
	}
	return freq.map(function (d) {return Math.round(d*100)/100;});
};


function make_gt_chart(gt){
	var tmp_data = [];
	var tmp_trace = ['x'];
	var tmp_colors = {};
	tmp_data.push(tmp_trace.concat(pivots));
	gt.forEach(function (d, i) {
		var region = d[0];
		var genotype = d[1];
		var freq = get_frequencies(region, genotype);
		console.log(region+' '+genotype);
		if (d3.max(freq)>0) {
			var tmp_trace = genotype.toString().replace(/,/g, ', ');
			if (region != "global") {
				tmp_trace = region + ':\t' + tmp_trace;
			}
			tmp_data.push([tmp_trace].concat(freq));
			tmp_colors[tmp_trace] = genotypeColors[i];
		}
	});
	console.log(tmp_colors);
	gt_chart.load({
       	columns: tmp_data,
       	unload: true
	});
	gt_chart.data.colors(tmp_colors);
}

function addClade(d) {
	if (typeof gt_chart != "undefined"){
		var plot_data = [['x'].concat(rootNode["pivots"])];
		var reg = "global";
		if ((typeof d.target.freq !="undefined" )&&(d.target.freq[reg] != "undefined")){
			plot_data[plot_data.length] = [reg].concat(d.target.freq[reg]);
		}
		if (plot_data.length > 1) {
			if (plot_data[1][0] == "global") {
				plot_data[1][0] = "clade";
			}
		}
		gt_chart.load({
	       	columns: plot_data
		});
	}
}

function removeClade() {
	if (typeof gt_chart != "undefined"){
		gt_chart.unload({
	       	ids: ["clade"]
		});
	}
}

width = parseInt(d3.select(".freqplot-container").style("width"), 10);
height = 250;
var position = "inset";

var gt_chart = c3.generate({
	bindto: '#gtchart',
	size: {width: width-10, height: height},
	onresize: function() {
		width = parseInt(d3.select(".freqplot-container").style("width"), 10);
		height = 250;
		gt_chart.resize({height: height, width: width});
	},
	legend: {
		position: position,
		inset: {
    		anchor: 'top-right',
    		x: 10,
    		y: -15,
    		step: 1
    	}
	},
  	color: {
        pattern: genotypeColors
    },
	axis: {
		y: {
			label: {
				text: 'frequency',
				position: 'outer-middle'
			},
			tick: {
				values: [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
				outer: false
			},
            min: 0,
			max: 1
		},
		x: {
			label: {
				text: 'time',
				position: 'outer-center'
			},
			tick: {
				values: time_ticks,
				outer: false
			}
		}
	},
	data: {
		x: 'x',
		columns: [],
	}
});

function contains(arr, obj) {
    for(var i=0; i<arr.length; i++) {
        if (arr[i] == obj) return true;
    }
}

d3.json(path + file_prefix + "frequencies.json", function(error, json){
	console.log(error);
	frequencies = json;
	pivots= frequencies["mutations"]["global"]["pivots"].map(function (d) {return Math.round(parseFloat(d)*100)/100;});
	var ticks = [Math.round(pivots[0])];
	var tick_step = Math.round((pivots[pivots.length-1]-pivots[0])/6*10)/10;
	while (ticks[ticks.length-1]<pivots[pivots.length-1]){
		ticks.push(Math.round((ticks[ticks.length-1]+tick_step)*10)/10);
	}
	//gt_chart.axis.x.values = ticks;
	/**
		parses a genotype string into region and positions
	**/

	var chart_data = {}
	var chart_types = {}
	var chart_xaxis = {}
	var posToAA = {};
	var ymin = 0;
	var xmax = 0;
	if (typeof genome_annotation !== 'undefined') {
		for (x in genome_annotation){
			chart_data['x'+x+'anno'] = genome_annotation[x][1];
			chart_data[x+'anno'] = genome_annotation[x][0].map(function(d) {return -0.1*d;});
			if (ymin>chart_data[x+'anno'][0]){
				ymin = chart_data[x+'anno'][0];
			}
			chart_types[x+'anno'] = 'line';
			chart_xaxis[x+'anno'] = 'x'+x+'anno';
		}
		ymin-=0.08;
	}

	for (gene in frequencies["entropy"]){
		chart_data[gene]=[];
		chart_data['x'+gene]=[];
		chart_types[gene]='bar';
		chart_xaxis[gene]='x'+gene;
		var offset = frequencies['location'][gene][0];
		for (var ii=0;ii<frequencies["entropy"][gene].length;ii+=1){
			if (Math.round(10000*frequencies["entropy"][gene][ii][1])/10000>0.05){
				chart_data[gene].push(Math.round(10000*frequencies["entropy"][gene][ii][1])/10000);
				chart_data['x'+gene].push(ii*3+offset);
				posToAA[ii*3+offset] = [gene, ii];
				if ((ii*3+offset)>xmax) {xmax = (ii*3+offset);}
			}
		}
	}
	var entropy_chart = c3.generate({
		bindto: '#entropy',
		size: {width: width-10, height: height},
		onresize: function() {
			width = parseInt(d3.select(".entropy-container").style("width"), 10);
			height = 250;
			entropy_chart.resize({height: height, width: width});
		},
		legend: {show: false},
		color: {pattern: ["#AAA"]},
		axis: {
			y: {
				label: {
					text: 'variability',
					position: 'outer-middle'
				},
				tick: {
					values: [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6],
					outer: false
				},
				min:ymin,
			},
			x: {
				label: {
					text: 'position',
					position: 'outer-center'
				},
				tick: {
					outer: false,
					values: ([1,2,3,4,5]).map(function (d){
						var dec = Math.pow(10,Math.floor(Math.log10(xmax/5)))
						var step = dec*Math.floor(xmax/5/dec);
						return d*step;
					})
				}
			},
		},
		data: {
			xs: chart_xaxis,
			json: chart_data,
			types: chart_types,
			onclick: function (d,i) {
            	gene = posToAA[d.x][0];
            	var pos = posToAA[d.x][1];
				if (frequencies["entropy"][gene][pos][2].length>1){
					var tmp = [];
					for (var ii=0;ii<frequencies["entropy"][gene][pos][2].length;ii+=1){
						tmp.push(["global",d.x+frequencies["entropy"][gene][pos][2][ii]]);
					}
					colorBy = "genotype";
					console.log("color by genotype: "+gene + ' ' + pos)
					colorByGenotypePosition([[gene, pos]]);
					d3.select("#gt-color").property("value", gene + ':' + (pos+1));
				}
		    },
		    onmouseover: function (d){
		    	document.body.style.cursor = "pointer";
		    },
		    onmouseout: function (d){
		    	document.body.style.cursor = "default";
		    },
			labels:{
				format:function (v, id, i, j){
					if ((typeof id !="undefined")&&(id.substring(id.length-4)=='anno')&&(i==1)){
						return id.substring(0,id.length-4);
					}else{return '';}
				}
			},
		},
		bar: {width: 2},
	    grid: {
    	    y: {
        	    lines: [{value: 0}]
        	},
        	focus:{
        		show:false
        	}
    	},
	    tooltip: {
	        format: {
	            title: function (d) {
	            	if (typeof posToAA[d] != "undefined"){
		            	var gene = posToAA[d][0];
		            	var pos = posToAA[d][1];
		            	return gene + ' codon ' + (pos+1) + frequencies["entropy"][gene][pos][2].join(",");
		            }else{ return d;}},
	            value: function (value, ratio, id) {
	                return id.substring(id.length-4)=='anno'?"start/stop":"Variability: "+value;
	            }
	        }
		},
	});

	d3.select("#plotfreq")
		.on("click", function (){
			gt = parse_gt_string(document.getElementById("gtspec").value);
			make_gt_chart(gt);
		});
	make_gt_chart(parse_gt_string(document.getElementById("gtspec").value));
});

function addSequence(current_seq, current_seq_name, seqs, all_names){
    if (current_seq_name==""){
        current_seq_name="input sequence";
    }
    var name_count = 0;
    for (var tmpii=0; tmpii<all_names; tmpii++){
        if (all_names[ii]==current_seq_name){
            name_count++;
        }
    }
    if (name_count){
        suffix=" "+(name_count+1);
    }else{suffix="";}
    all_names.push(current_seq_name);
    seqs[current_seq_name+suffix]=current_seq;    
}

function parseSequences(){
    var lines = document.getElementById('seqinput').value.split('\n');
    var seqs = {};
    var unmatched = [];
    var closest_nodes = {};
    var current_seq_name = "";
    var current_seq = "";
    var seq_names = [];
    var suffix;
    for (var li=0; li<lines.length; li++){
        if (lines[li][0]=='>'){
            if (current_seq.length){
                addSequence(current_seq, current_seq_name, seqs, seq_names);  
            }
            current_seq_name = lines[li].substring(1,lines[li].length);
            current_seq = "";
        }else{
            current_seq += lines[li].toUpperCase().replace(/[^ACGTWRN]/g,"");
        }
    }
    if (current_seq.length){
        addSequence(current_seq, current_seq_name, seqs, seq_names);  
    }
    for (current_seq_name in seqs){
        var tmpclade = locateSequence(current_seq_name, seqs[current_seq_name]);
        if (tmpclade!=null){
            if (typeof closest_nodes[tmpclade]=="undefined"){closest_nodes[tmpclade]=[current_seq_name];}
            else{closest_nodes[tmpclade].push(current_seq_name);}
        }else{
            unmatched.push(current_seq_name);
        }
    }
    markInTreeSeqSearch(closest_nodes);
    if (unmatched.length){
        console.log(unmatched);
        var tmp_str = "";
        for (var ii=0; ii<unmatched.length; ii++){ tmp_str+=unmatched[ii].substring(0,30)+"\n";}
        window.alert("No close match was found for \n" + tmp_str + "\nMaybe gapped or not recent isolate from current lineage?");
    }
}

function locateSequence(name, seq){
    var mutations, olap_start, olap_end;
    console.log('Provided sequence: '+ name +': ' + seq.substring(0,20)+'....');
    tmp = alignToRoot(seq);
    olap_start=tmp[0]; olap_end=tmp[1]; mutations=tmp[2];
    if (olap_start==null){
        return null;
    }else{
        console.log("start, end:", olap_start, olap_end);
        console.log("mutations:", mutations);
        var bestClade = findClosestClade(mutations);
        return bestClade;
    }
}

function findClosestClade(mutations){
    var bestClade=-1, bestScore=0;
    var tmpScore=0;
    var searchClades = tips.map(function(d){return d.clade;});

    for (ci=0; ci<searchClades.length; ci++){
        clade = searchClades[ci];
        tmpScore=0;
        for (mut in mutations){
            if (stateAtPosition(clade, 'nuc', mut)==mutations[mut]){
                tmpScore++;
            }
        }
        if (clade!="root") {
            tmpScore -= 0.5*Object.keys(cladeToSeq[clade]['nuc']).length;
        }
        if (tmpScore>bestScore){
            bestScore=tmpScore;
            bestClade=clade;
        }
    }
    console.log("best match:",bestClade);
    return bestClade;
}


function alignToRoot(seq){
    var rootSeq = cladeToSeq["root"]["nuc"];
    var shift = 0;
    var max_score = 0.0, max_shift;

    for(shift=0; shift<seq.length-30;shift++){
        var tmp_score = 0;
        var olaplen=Math.min(seq.length-shift, rootSeq.length);
        for (var pos=0; pos<olaplen; pos++){
            if (seq[pos+shift]==rootSeq[pos]){
                tmp_score++;
            }
        }
        tmp_score*=1.0/olaplen;
        if (tmp_score>max_score){
            max_score=tmp_score;
            max_shift=-shift;
        }
    }

    for(shift=0; shift<rootSeq.length-30;shift++){
        var tmp_score = 0;
        var olaplen=Math.min(rootSeq.length-shift, seq.length);
        for (var pos=0; pos<olaplen; pos++){
            if (seq[pos]==rootSeq[shift+pos]){
                tmp_score++;
            }
        }
        tmp_score*=1.0/olaplen;
        if (tmp_score>max_score){
            max_score=tmp_score;
            max_shift=shift;
        }
    }
    console.log("best shift: ",max_shift, " score: ",max_score);
    if (max_score>0.9){
        var mutations = {};
        if (max_shift<0){
            var olaplen=Math.min(seq.length+max_shift, rootSeq.length);
            var olap_start = 0;
            var olap_end = olaplen;
        }else{
            var olaplen=Math.min(rootSeq.length-max_shift, seq.length);
            var olap_start = max_shift;
            var olap_end = max_shift+olaplen;
        }
        for (var pos=olap_start; pos<olap_end; pos++){
            if (rootSeq[pos]!=seq[pos-max_shift]){
                mutations[pos]=seq[pos-max_shift];
            }
        }
        return [olap_start, olap_end, mutations];
    }else{
        console.log("no good match");
        return [null, null, null];
    }
}

// highlight clades in tree
function markInTreeSeqSearch(clades){
    var userSeqs = nodes.filter(function(d){
        var tmp=0; for (var clade in clades){tmp+= (d.clade==clade);} return tmp>0;});

    for (var mi=0; mi<userSeqs.length; mi++){
        userSeqs[mi].matches = clades[userSeqs[mi].clade];
    }

    treeplot.selectAll('.seqmatch').data(userSeqs)
        .enter()
        .append('text')
        .attr("class", "seqmatch")
        .text(function(d) { console.log(d.strain); return '\uf069'; })
        .on('mouseover', function(d) {
            matchTooltip.show(d, this);
        })
        .on('mouseout', matchTooltip.hide);
    styleHighlight();
}

function markInTreeStrainSearch(tip){
    treeplot.selectAll('.strainmatch').data([tip])
        .enter()
        .append('text')
        .attr("class", "strainmatch")        
        .text(function(d) { console.log(d.strain); return '\uf069'; })
        .on('mouseover', function(d) {
            virusTooltip.show(d, this);
        })
        .on('mouseout', virusTooltip.hide);
    styleHighlight();
}

function styleHighlight(){
    treeplot.selectAll('.seqmatch')
        .attr('text-anchor', 'middle')
        .attr('dominant-baseline', 'central')
        .style("font-size", "24px")
        .style('font-family', 'FontAwesome')
        .style("fill", "#555555")
        .attr("x", function(d) { return d.x; })
        .attr("y", function(d) { return d.y; })
        .style("cursor", "default"); 
    treeplot.selectAll('.strainmatch')
        .attr('text-anchor', 'middle')
        .attr('dominant-baseline', 'central')
        .style("font-size", "24px")
        .style('font-family', 'FontAwesome')
        .style("fill", "#555555")
        .attr("x", function(d) { return d.x; })
        .attr("y", function(d) { return d.y; })
        .style("cursor", "default");              
}

// callback to highlight the result of a search by strain name
var searchEvent;
function highlightStrainSearch(tip) {
    var strainName = (tip.strain).replace(/\//g, "");
    d3.select("#"+strainName)
        .call(function(d) {
            markInTreeStrainSearch(tip);
            virusTooltip.show(tip, d[0][0]);
        });
}

var strainSearchEvent;
d3.select('#seqinput').on('keyup', function(){
        if (typeof strainSearchEvent != "undefined"){clearTimeout(strainSearchEvent);}
        strainSearchEvent = setTimeout(parseSequences, 100);
    }); 

d3.select('#searchinputclear').on('click', function (){
    treeplot.selectAll('.seqmatch').data([]).exit().remove();
    treeplot.selectAll('.strainmatch').data([]).exit().remove();    
    document.getElementById('seqinput').value = "";
    document.getElementById('bp-input').value = "";
	virusTooltip.hide();
    });

